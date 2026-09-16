#!/usr/bin/env python3


import sys
import tempfile
from pathlib import Path
from shutil import which
from subprocess import PIPE, Popen

from yaml import safe_load

_TIBERIUS_REPO = "https://github.com/Gaius-Augustus/Tiberius"
_YAML_KEYS = ["ncbi_tax_id", "softmasking"]
_GIT = Path(which("git"))


def get_model_dict(yaml_file):
    with open(yaml_file, "rt") as f:
        model_config = safe_load(f)
        model_dict = {k: model_config.get(k) for k in _YAML_KEYS}
        model_dict["model_cfg"] = yaml_file.name
    return model_dict


def checkout_tag(tool_repo: Path, tool_version: str) -> bool:
    # pull the tag
    git_cmd = [
        _GIT,
        "-C",
        tool_repo,
        "checkout",
        tool_version,
    ]
    with Popen(git_cmd, stdout=PIPE) as proc:
        print(proc.stdout.read(), file=sys.stderr)

    return True


def get_release(tool_repo: Path) -> str:
    git_cmd = [
        _GIT,
        "-C",
        tool_repo,
        "describe",
        "--tags",
        "--abbrev=0",
        "--match=v*",
        "--exclude=v*-",
    ]

    with Popen(git_cmd, stdout=PIPE, stderr=PIPE) as proc:
        out, err = proc.communicate()

    print(err.decode(), file=sys.stderr)

    return out.decode().rstrip()


def pull_repo() -> Path:

    tempdir = tempfile.mkdtemp()

    # pull the tag
    pull_command = [
        _GIT,
        "clone",
        _TIBERIUS_REPO,
        tempdir,
    ]
    with Popen(pull_command, stdout=PIPE) as proc:
        print(proc.stdout.read(), file=sys.stderr)

    return Path(tempdir)


def main():
    tool_repo = pull_repo()
    tiberius_release = get_release(tool_repo)
    print(f"Getting models for tiberius_release {tiberius_release}", file=sys.stderr)
    _ = checkout_tag(tool_repo, tiberius_release)

    yaml_files = Path(tool_repo, "model_cfg").glob("*.yaml")
    model_dicts = [get_model_dict(x) for x in yaml_files]

    print(
        f"# Tiberius tag {tiberius_release} non-softmasking models from {_TIBERIUS_REPO}",
        file=sys.stdout,
    )
    for model_dict in model_dicts:
        model_cfg = model_dict.get("model_cfg")
        if model_dict.get("softmasking") == False:
            ncbi_tax_id = model_dict.get("ncbi_tax_id")
            print(f"{ncbi_tax_id}\t{model_cfg}", file=sys.stdout)
        else:
            print(f"Skipping softmasking model {model_cfg}", file=sys.stderr)


if __name__ == "__main__":

    main()
