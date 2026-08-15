#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pure-Python remote tabix point-query for bgzipped, tabix-indexed GWAS sumstats
served over HTTP with byte-range support (e.g. FinnGen R12 on Google Storage).

No pysam / C compiler required. Uses only urllib + gzip + struct.

Strategy for a POINT query (chrom, pos):
  * parse the .tbi index (small; downloaded/cached once per endpoint)
  * use the LINEAR index (16 kb windows) to get the smallest virtual file
    offset of records overlapping the query window
  * HTTP byte-range GET the BGZF blocks from that offset, gunzip forward
    (BGZF is concatenated standard gzip members), and scan TSV rows until
    the record position exceeds the query -> return matching rows.

Only the bytes around each queried SNP are transferred, so all ~2,469 FinnGen
endpoints can be scanned without downloading the multi-hundred-MB files.
"""
import io
import gzip
import zlib
import struct
import urllib.request


def bgzf_decompress(blob, max_out=8_000_000):
    """Decode concatenated BGZF blocks from `blob`, tolerating a truncated
    final block. Returns (decompressed_bytes, consumed_compressed_bytes)."""
    out = bytearray()
    off = 0
    n = len(blob)
    while off + 18 <= n and len(out) < max_out:
        if blob[off:off + 2] != b"\x1f\x8b":
            break
        xlen = struct.unpack_from("<H", blob, off + 10)[0]
        # find BSIZE in extra field (subfield SI1='B',SI2='C',SLEN=2)
        bsize = None
        ext = off + 12
        end = ext + xlen
        p = ext
        while p + 4 <= end:
            si1, si2, slen = blob[p], blob[p + 1], struct.unpack_from("<H", blob, p + 2)[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack_from("<H", blob, p + 4)[0]
                break
            p += 4 + slen
        if bsize is None:
            break
        block_len = bsize + 1
        if off + block_len > n:
            break  # truncated final block -> stop cleanly
        cdata_start = off + 12 + xlen
        cdata_end = off + block_len - 8
        try:
            out += zlib.decompress(blob[cdata_start:cdata_end], -15)
        except zlib.error:
            break
        off += block_len
    return bytes(out), off


def _http_range(url, start, length, timeout=60, retries=3):
    end = start + length - 1
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    last = None
    for _ in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
    raise last


def _http_get(url, timeout=120, retries=3):
    last = None
    for _ in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
    raise last


class TabixIndex:
    """Parsed .tbi index: chrom-name -> tid and per-tid 16kb linear index."""

    def __init__(self, tbi_bytes):
        raw = gzip.decompress(tbi_bytes)
        self._parse(raw)

    def _parse(self, d):
        off = 0

        def u32():
            nonlocal off
            v = struct.unpack_from("<i", d, off)[0]
            off += 4
            return v

        def u64():
            nonlocal off
            v = struct.unpack_from("<q", d, off)[0]
            off += 8
            return v

        magic = d[off:off + 4]
        off += 4
        assert magic == b"TBI\x01", f"bad tbi magic {magic!r}"
        n_ref = u32()
        self.fmt = u32()
        self.col_seq = u32()
        self.col_beg = u32()
        self.col_end = u32()
        self.meta = u32()
        self.skip = u32()
        l_nm = u32()
        names_blob = d[off:off + l_nm]
        off += l_nm
        names = names_blob.split(b"\x00")
        names = [n.decode() for n in names if n]
        self.name2tid = {n: i for i, n in enumerate(names)}

        self.linear = []  # per tid: list of virtual offsets (uint64)
        for _ in range(n_ref):
            n_bin = u32()
            for _b in range(n_bin):
                _bin_id = struct.unpack_from("<I", d, off)[0]
                off += 4
                n_chunk = u32()
                off += 16 * n_chunk  # skip chunk begin/end pairs
            n_intv = u32()
            intv = list(struct.unpack_from("<%dq" % n_intv, d, off))
            off += 8 * n_intv
            self.linear.append(intv)

    def tid(self, chrom):
        if chrom in self.name2tid:
            return self.name2tid[chrom]
        # tolerate chr-prefix mismatch
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        return self.name2tid.get(alt)

    def voffset(self, chrom, pos):
        """virtual offset to start scanning for 1-based position `pos`."""
        tid = self.tid(chrom)
        if tid is None or tid >= len(self.linear):
            return None
        intv = self.linear[tid]
        if not intv:
            return None
        w = (pos - 1) >> 14  # 16 kb window, 0-based coordinate
        if w >= len(intv):
            w = len(intv) - 1
        # first non-zero offset at or before window (0 means "no record yet")
        vo = intv[w]
        while vo == 0 and w > 0:
            w -= 1
            vo = intv[w]
        return vo if vo else intv[0]


class RemoteTabix:
    def __init__(self, gz_url, tbi_bytes=None):
        self.url = gz_url
        if tbi_bytes is None:
            tbi_bytes = _http_get(gz_url + ".tbi")
        self.idx = TabixIndex(tbi_bytes)
        self.col_seq = self.idx.col_seq
        self.col_beg = self.idx.col_beg

    def query(self, chrom, pos, fetch=262144, max_scan=4_000_000):
        """Return list of TSV field-lists at exactly (chrom, pos). 1-based pos."""
        vo = self.idx.voffset(chrom, pos)
        if vo is None:
            return []
        coffset = vo >> 16
        uoffset = vo & 0xFFFF
        blob = _http_range(self.url, coffset, fetch)
        if not blob:
            return []
        # BGZF == concatenated gzip blocks; decode forward from a block start
        data, _ = bgzf_decompress(blob, max_out=max_scan)
        stream = io.BytesIO(data[uoffset:])
        seqi = self.col_seq - 1
        begi = self.col_beg - 1
        out = []
        want_seq = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        need_more = False
        for line in stream:
            try:
                s = line.decode()
            except UnicodeDecodeError:
                continue
            if not s.endswith("\n"):
                need_more = True  # possibly truncated last line
                break
            if s[0] in "#":
                continue
            f = s.rstrip("\n").split("\t")
            if len(f) <= begi:
                continue
            if f[seqi] not in want_seq:
                continue
            try:
                p = int(f[begi])
            except ValueError:
                continue
            if p == pos:
                out.append(f)
            elif p > pos:
                need_more = False
                break
        # if we ran out of fetched bytes before reaching pos, widen once
        if need_more and fetch < 4_000_000:
            return self.query(chrom, pos, fetch=fetch * 4, max_scan=max_scan)
        return out

    def region(self, chrom, start, end, chunk=1_048_576, max_bytes=64_000_000):
        """Return all TSV field-lists with start <= pos <= end (1-based).

        Fetches BGZF blocks forward from the linear-index offset for `start`,
        transparently continuing to pull more compressed bytes until a record
        with pos > end is seen (or max_bytes transferred)."""
        vo = self.idx.voffset(chrom, start)
        if vo is None:
            return []
        coffset = vo >> 16
        uoffset = vo & 0xFFFF
        seqi = self.col_seq - 1
        begi = self.col_beg - 1
        want_seq = {chrom, chrom[3:] if chrom.startswith("chr") else "chr" + chrom}
        out = []
        pending = b""          # carry an unterminated trailing line between fetches
        fetched = 0
        first = True
        while fetched < max_bytes:
            blob = _http_range(self.url, coffset + fetched, chunk)
            if not blob:
                break
            fetched += len(blob)
            data, consumed = bgzf_decompress(blob, max_out=16_000_000)
            if not data:
                break
            if first:
                data = data[uoffset:]
                first = False
            buf = pending + data
            nl = buf.rfind(b"\n")
            if nl < 0:
                pending = buf
                continue
            pending = buf[nl + 1:]
            done = False
            for line in buf[:nl].split(b"\n"):
                if not line or line[:1] == b"#":
                    continue
                try:
                    f = line.decode().split("\t")
                except UnicodeDecodeError:
                    continue
                if len(f) <= begi or f[seqi] not in want_seq:
                    continue
                try:
                    p = int(f[begi])
                except ValueError:
                    continue
                if p < start:
                    continue
                if p > end:
                    done = True
                    break
                out.append(f)
            if done:
                break
        return out


def _partial_gunzip(blob, max_scan):
    out = bytearray()
    bio = io.BytesIO(blob)
    while len(out) < max_scan:
        try:
            g = gzip.GzipFile(fileobj=bio)
            chunk = g.read()
        except (OSError, EOFError):
            break
        if not chunk:
            break
        out += chunk
        # advance bio past this member is handled by GzipFile; if stuck, stop
        if bio.tell() >= len(blob):
            break
    return bytes(out)
