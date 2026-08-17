#!/usr/bin/env python
"""
auto_sync.py — keep the published Human Plasma Immune Atlas in step with this project.

One command rebuilds the app data from the current analysis output and pushes
*both* repositories:

    Huggeface/  ->  Hugging Face Space  jianlizhao/Human-Plasma-Immune-Atlas
                ->  GitHub              jimmyuab/Human-Plasma-Immune-Atlas
    github/     ->  GitHub              jimmyuab/Human-Plasma-Immune-Atlas-Code

Usage
-----
    python auto_sync.py --save-hf-token       # store a Hugging Face WRITE token (once)
    python auto_sync.py --save-github-token   # store a GitHub PAT with `repo` scope (once)
    python auto_sync.py --check               # who am I logged in as? what is publishable?
    python auto_sync.py                       # rebuild + commit + push everything
    python auto_sync.py --no-build            # skip the data rebuild
    python auto_sync.py --install-schedule    # run it automatically every day at 20:00
    python auto_sync.py --remove-schedule

Where the login is kept
-----------------------
Hugging Face : .secrets/hf_token.txt      (git-ignored)  + the huggingface_hub
                                          credential store, so `hf`/`huggingface-cli`
                                          and the Python client stay logged in too.
GitHub       : .secrets/github_token.txt  (git-ignored)
Either can be overridden for a single run with the HF_TOKEN / GITHUB_TOKEN env vars.

Tokens belong to you. Create them yourself at
    https://huggingface.co/settings/tokens   (role: write)
    https://github.com/settings/tokens       (scope: repo)
and paste them in when --save-*-token asks. Nothing in this project ever
generates, guesses or transmits a credential anywhere except to the service it
authenticates.
"""
from __future__ import annotations

import argparse
import datetime
import getpass
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.abspath(__file__))
SECRETS = os.path.join(ROOT, ".secrets")
APP = os.path.join(ROOT, "Huggeface")
CODE = os.path.join(ROOT, "github")

HF_USER, HF_SPACE = "jianlizhao", "Human-Plasma-Immune-Atlas"
GH_USER = "jimmyuab"
GH_APP_REPO = "Human-Plasma-Immune-Atlas"
GH_CODE_REPO = "Human-Plasma-Immune-Atlas-Code"

HF_URL = f"https://huggingface.co/spaces/{HF_USER}/{HF_SPACE}"
GH_APP_URL = f"https://github.com/{GH_USER}/{GH_APP_REPO}.git"
GH_CODE_URL = f"https://github.com/{GH_USER}/{GH_CODE_REPO}.git"

TASK_NAME = "PlasmaImmuneAtlas-AutoSync"


# ---------------------------------------------------------------- plumbing --
def git_env(repo: str) -> dict:
    """The I: drive does not record ownership, so git refuses with 'dubious
    ownership'. Scope the exception to this process instead of writing it into
    the machine's global git config."""
    return dict(os.environ,
                GIT_CONFIG_COUNT="1",
                GIT_CONFIG_KEY_0="safe.directory",
                GIT_CONFIG_VALUE_0=repo.replace("\\", "/"))


def redact(text: str) -> str:
    for tok in ("hf_", "ghp_", "github_pat_"):
        if tok in text:
            head, _, _ = text.partition(tok)
            text = head + tok + "***"
    return text


def sh(cmd, repo: str, check=True, quiet=False):
    if not quiet:
        print("  $ " + redact(cmd if isinstance(cmd, str) else " ".join(cmd)))
    return subprocess.run(cmd, cwd=repo, shell=isinstance(cmd, str), env=git_env(repo),
                          check=check, text=True, capture_output=True)


def step(msg: str):
    print(f"\n=== {msg} ===")


# ------------------------------------------------------------ credentials --
def read_token(env_name: str, filename: str) -> str | None:
    v = os.environ.get(env_name)
    if v and v.strip():
        return v.strip()
    p = os.path.join(SECRETS, filename)
    if os.path.exists(p):
        v = open(p).read().strip()
        if v:
            return v
    return None


def hf_token() -> str | None:
    t = read_token("HF_TOKEN", "hf_token.txt")
    if t:
        return t
    try:                                    # fall back to an existing hf login
        from huggingface_hub import get_token
        return get_token()
    except Exception:
        return None


def github_token() -> str | None:
    return read_token("GITHUB_TOKEN", "github_token.txt")


def save_hf_token():
    tok = getpass.getpass("Paste your Hugging Face WRITE token (input hidden): ").strip()
    if not tok:
        sys.exit("nothing entered — aborted")
    try:
        from huggingface_hub import whoami
        who = whoami(tok)
    except Exception as e:
        sys.exit(f"Hugging Face rejected that token: {e}")
    name = who.get("name")
    if name != HF_USER:
        print(f"  ! the token belongs to '{name}', not '{HF_USER}' — pushes to the "
              f"{HF_USER} Space will fail unless you have write access to it.")
    os.makedirs(SECRETS, exist_ok=True)
    path = os.path.join(SECRETS, "hf_token.txt")
    with open(path, "w") as fh:
        fh.write(tok + "\n")
    try:                                    # also keep the CLI / client logged in
        from huggingface_hub import login
        login(token=tok, add_to_git_credential=False)
        print("  huggingface_hub credential store: logged in")
    except Exception as e:
        print(f"  (could not write the huggingface_hub store: {e})")
    print(f"  saved -> {path}   (git-ignored)")
    print(f"  logged in as: {name}")


def save_github_token():
    tok = getpass.getpass("Paste your GitHub PAT with 'repo' scope (input hidden): ").strip()
    if not tok:
        sys.exit("nothing entered — aborted")
    os.makedirs(SECRETS, exist_ok=True)
    path = os.path.join(SECRETS, "github_token.txt")
    with open(path, "w") as fh:
        fh.write(tok + "\n")
    print(f"  saved -> {path}   (git-ignored)")


# ------------------------------------------------------------------ checks --
def check():
    step("credentials")
    ht = hf_token()
    if ht:
        try:
            from huggingface_hub import whoami
            print(f"  Hugging Face : logged in as {whoami(ht).get('name')}")
        except Exception as e:
            print(f"  Hugging Face : token present but REJECTED — {e}")
    else:
        print("  Hugging Face : no token  ->  python auto_sync.py --save-hf-token")
    if github_token():
        print("  GitHub       : personal access token present")
    else:
        r = subprocess.run(["git", "credential", "fill"], input="protocol=https\nhost=github.com\n\n",
                           text=True, capture_output=True, env=dict(os.environ, GIT_TERMINAL_PROMPT="0"))
        print("  GitHub       : " + ("Credential Manager has a saved github.com login"
                                     if "password=" in r.stdout else
                                     "not logged in  ->  the first push opens a browser sign-in"))

    step("the Space")
    if ht:
        try:
            from huggingface_hub import HfApi
            info = HfApi().space_info(f"{HF_USER}/{HF_SPACE}", token=ht)
            print(f"  {HF_USER}/{HF_SPACE}: exists, visibility="
                  f"{'PRIVATE' if info.private else 'PUBLIC'}, sdk={info.sdk}")
            if info.private:
                print("  ! The Space is PRIVATE. The atlas is meant to be usable with no "
                      "login, so flip it to Public:\n"
                      f"    {HF_URL}/settings  ->  Change visibility  ->  Public")
        except Exception as e:
            print(f"  {HF_USER}/{HF_SPACE}: not reachable — {e}")
            print(f"  Create it: https://huggingface.co/new-space  ->  name {HF_SPACE}, "
                  f"SDK Gradio, visibility PUBLIC")

    step("local repositories")
    for label, repo in (("app  (Huggeface/)", APP), ("code (github/)", CODE)):
        if not os.path.isdir(os.path.join(repo, ".git")):
            print(f"  {label}: not a git repository")
            continue
        dirty = sh("git status --porcelain", repo, check=False, quiet=True).stdout.strip()
        head = sh(["git", "log", "-1", "--format=%h %s"], repo, check=False, quiet=True).stdout.strip()
        remotes = sh("git remote", repo, check=False, quiet=True).stdout.split()
        print(f"  {label}: {head}")
        print(f"       uncommitted files: {len(dirty.splitlines()) if dirty else 0}"
              f" | remotes: {', '.join(remotes) if remotes else 'none'}")


# ------------------------------------------------------------------- push ---
def ensure_identity(repo: str):
    if not sh("git config user.email", repo, check=False, quiet=True).stdout.strip():
        sh('git config user.email "atlas@plasma-immunome.local"', repo, quiet=True)
        sh('git config user.name "Human Plasma Immune Atlas"', repo, quiet=True)
    # Git Credential Manager remembers the GitHub browser login, so no personal
    # access token is needed for GitHub. Set it per-repo, never globally.
    sh("git config --local credential.helper manager", repo, check=False, quiet=True)
    sh("git config --local credential.https://github.com.provider github", repo,
       check=False, quiet=True)


def commit(repo: str, message: str) -> bool:
    ensure_identity(repo)
    sh("git add -A", repo, quiet=True)
    if not sh("git diff --cached --quiet", repo, check=False, quiet=True).returncode:
        print("  nothing new to commit")
        return False
    sh(["git", "commit", "-m", message], repo)
    return True


def gh_auth(url: str, tok: str | None) -> str:
    """A PAT if one was saved, otherwise the plain URL so Git Credential Manager
    supplies the browser login it already remembers."""
    return url.replace("https://", f"https://{GH_USER}:{tok}@") if tok else url


def existing_remote(repo: str, name: str) -> str | None:
    """Whatever the repo is already wired to wins, so renaming the GitHub repo
    does not require editing this script."""
    r = sh(f"git remote get-url {name}", repo, check=False, quiet=True)
    url = r.stdout.strip()
    return url if url.startswith("https://") else None


def push(repo: str, name: str, auth_url: str, plain_url: str) -> bool:
    if sh(f"git remote get-url {name}", repo, check=False, quiet=True).returncode:
        sh(["git", "remote", "add", name, auth_url], repo, quiet=True)
    else:
        sh(["git", "remote", "set-url", name, auth_url], repo, quiet=True)
    r = sh(f"git push -u {name} HEAD:main", repo, check=False)
    rejected = any(k in (r.stderr or "") + (r.stdout or "")
                   for k in ("non-fast-forward", "fetch first", "Updates were rejected"))
    if r.returncode and name == "hf" and rejected:
        # A freshly created Space is seeded with its own README commit, so the first
        # push diverges. The Space is a pure mirror of this repo — nothing is authored
        # there — so overwriting that seed commit loses no work.
        print("  the Space has its own initial commit; overwriting it with this repo")
        r = sh(f"git push -u --force {name} HEAD:main", repo, check=False)
    sh(["git", "remote", "set-url", name, plain_url], repo, quiet=True)  # never persist the token
    if r.returncode:
        print(redact(r.stdout or ""))
        print(redact(r.stderr or ""))
        return False
    print(f"  pushed -> {plain_url}")
    return True


# ------------------------------------------------------------ code staging ---
# What goes into the public code repository, as (source under ROOT, destination
# under github/). Everything is copied newest-wins; nothing is deleted, so an old
# result that is no longer produced stays in the repo rather than vanishing silently.
CODE_MAP = [
    ("src",                                    "src"),
    ("09_tables",                              "results/tables"),
    ("08_figures",                             "figures"),
    ("10_manuscript",                          "manuscript"),
]
# individual result tables are named explicitly: 06_genetic_causality also holds the
# 615 MB cis_MR_ALL scan and the per-disease panphenome/ tree, neither of which can go
# to GitHub (100 MB file limit) and both of which src/41 regenerates.
CODE_RESULTS = os.path.join(ROOT, "06_genetic_causality")
CODE_SKIP_DIRS = {"panphenome", "_ppt_cache", "_ppt_render", "__pycache__"}
# dated working copies of files the repo already carries under their canonical name
CODE_SKIP_FILES = {"0815Methodology.docx"}
CODE_MAX_MB = 90


def stage_code_repo() -> int:
    """Copy the current pipeline output into github/ before it is committed.

    Without this the daily scheduled run pushed whatever happened to have been
    copied across by hand, so a re-run of the analysis never reached GitHub — the
    repo silently froze while the project moved on.
    """
    import shutil
    copied = 0

    def newer(src, dst):
        return (not os.path.exists(dst)
                or os.path.getmtime(src) > os.path.getmtime(dst) + 1
                or os.path.getsize(src) != os.path.getsize(dst))

    def put(src, dst):
        nonlocal copied
        if os.path.getsize(src) > CODE_MAX_MB * 1024 * 1024:
            return                                    # over GitHub's file limit
        if not newer(src, dst):
            return
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        copied += 1

    for rel_src, rel_dst in CODE_MAP:
        base = os.path.join(ROOT, rel_src)
        if not os.path.isdir(base):
            continue
        for dirpath, dirnames, filenames in os.walk(base):
            dirnames[:] = [d for d in dirnames if d not in CODE_SKIP_DIRS]
            for f in filenames:
                if f.startswith("~$") or f.endswith(".pyc") or f in CODE_SKIP_FILES:
                    continue
                s = os.path.join(dirpath, f)
                put(s, os.path.join(CODE, rel_dst, os.path.relpath(s, base)))

    put(os.path.join(ROOT, "auto_sync.py"), os.path.join(CODE, "auto_sync.py"))

    for f in sorted(os.listdir(CODE_RESULTS)):
        s = os.path.join(CODE_RESULTS, f)
        if os.path.isfile(s) and f.endswith((".tsv", ".csv", ".md")):
            put(s, os.path.join(CODE, "results", "genetic_causality", f))

    print(f"  staged {copied} changed file(s) into github/")
    return copied


def sync(build: bool, message: str) -> bool:
    ok = True

    if build:
        step("rebuilding Huggeface/data from the current analysis output")
        if subprocess.run([sys.executable, os.path.join(APP, "build_data.py")], cwd=APP).returncode:
            sys.exit("data rebuild failed — fix that before publishing")

    ht, gt = hf_token(), github_token()

    step("app repository (Huggeface/)")
    commit(APP, message)
    if ht:
        ok &= push(APP, "hf", f"https://{HF_USER}:{ht}@huggingface.co/spaces/{HF_USER}/{HF_SPACE}",
                   HF_URL)
    else:
        print("  ! no Hugging Face token — skipped. python auto_sync.py --save-hf-token")
        ok = False
    app_url = existing_remote(APP, "github") or GH_APP_URL
    ok &= push(APP, "github", gh_auth(app_url, gt), app_url)

    step("code repository (github/)")
    stage_code_repo()
    commit(CODE, message)
    code_url = existing_remote(CODE, "origin") or GH_CODE_URL
    ok &= push(CODE, "origin", gh_auth(code_url, gt), code_url)

    print("\n" + "=" * 66)
    print("DONE" if ok else "FINISHED WITH WARNINGS — read the messages above")
    print(f"  Space : {HF_URL}")
    print(f"  App   : {GH_APP_URL[:-4]}")
    print(f"  Code  : {GH_CODE_URL[:-4]}")
    print("  The Space rebuilds itself ~2-4 min after the push.")
    print("=" * 66)
    return ok


# --------------------------------------------------------------- schedule ---
def install_schedule(at: str):
    cmd = f'"{sys.executable}" "{os.path.join(ROOT, "auto_sync.py")}"'
    r = subprocess.run(["schtasks", "/Create", "/TN", TASK_NAME, "/TR", cmd,
                        "/SC", "DAILY", "/ST", at, "/F"],
                       text=True, capture_output=True)
    print(r.stdout or r.stderr)
    if r.returncode:
        sys.exit("could not register the scheduled task")
    print(f"  '{TASK_NAME}' will rebuild and publish every day at {at}.")
    print(f"  Remove it with: python auto_sync.py --remove-schedule")


def remove_schedule():
    r = subprocess.run(["schtasks", "/Delete", "/TN", TASK_NAME, "/F"],
                       text=True, capture_output=True)
    print(r.stdout or r.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--save-hf-token", action="store_true")
    ap.add_argument("--save-github-token", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--message", "-m", default=None)
    ap.add_argument("--install-schedule", action="store_true")
    ap.add_argument("--at", default="20:00", help="time of day for --install-schedule")
    ap.add_argument("--remove-schedule", action="store_true")
    a = ap.parse_args()

    if a.save_hf_token:
        return save_hf_token()
    if a.save_github_token:
        return save_github_token()
    if a.install_schedule:
        return install_schedule(a.at)
    if a.remove_schedule:
        return remove_schedule()
    if a.check:
        return check()

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    sync(build=not a.no_build,
         message=a.message or f"Update Human Plasma Immune Atlas ({stamp})")


if __name__ == "__main__":
    main()
