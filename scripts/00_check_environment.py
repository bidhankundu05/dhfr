"""
Module 0: Environment check.

Verifies that the required Python packages import correctly and that the
external command-line tools used later in the pipeline (MAFFT, AutoDock
Vina, Open Babel, FastTree) are installed and on PATH. Prints a PASS/FAIL
summary table and, for any missing CLI tool, an OS-appropriate install
command.
"""

import platform
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# 1. Python package import checks
# ---------------------------------------------------------------------------

PACKAGE_CHECKS = [
    ("Biopython", "Bio", lambda mod: getattr(mod, "__version__", "unknown")),
    ("GEOparse", "GEOparse", lambda mod: getattr(mod, "__version__", "unknown")),
    ("pandas", "pandas", lambda mod: getattr(mod, "__version__", "unknown")),
    ("pubchempy", "pubchempy", lambda mod: getattr(mod, "__version__", "unknown")),
]


def check_packages():
    results = []
    for label, module_name, version_fn in PACKAGE_CHECKS:
        try:
            mod = __import__(module_name)
            version = version_fn(mod)
            results.append((label, "PASS", f"imported OK (version {version})"))
        except Exception as exc:  # noqa: BLE001 - want to report any import failure
            results.append((label, "FAIL", f"import failed: {exc}"))
    return results


# ---------------------------------------------------------------------------
# 2. Command-line tool checks
# ---------------------------------------------------------------------------

# (display name, executable name, version flag, fallback executable/version-flag or None)
# FastTree itself has no macOS Homebrew formula; VeryFastTree is an actively
# maintained, CLI-compatible reimplementation (same flags: -nt, -gtr, -quiet,
# -log, -pseudo, ...) and is accepted here as an explicit, disclosed substitute
# -- same spirit as the PDB/AlphaFold substitution for Module 5.
CLI_TOOLS = [
    ("MAFFT", "mafft", "--version", None),
    ("AutoDock Vina", "vina", "--version", None),
    ("Open Babel", "obabel", "-V", None),
    ("FastTree", "fasttree", "", ("veryfasttree", "--help")),
]

# Install commands per OS. FastTree has no --version flag and often prints
# its usage/version banner to stderr on a bare invocation.
INSTALL_COMMANDS = {
    "Darwin": {
        "mafft": "brew install mafft",
        "vina": "brew install autodock-vina",
        "obabel": "brew install open-babel",
        "fasttree": "brew install fasttree",
    },
    "Linux": {
        "mafft": "sudo apt-get install -y mafft",
        "vina": "sudo apt-get install -y autodock-vina  (or download from https://vina.scripps.edu)",
        "obabel": "sudo apt-get install -y openbabel",
        "fasttree": "sudo apt-get install -y fasttree",
    },
    "Windows": {
        "mafft": "conda install -c bioconda mafft",
        "vina": "conda install -c bioconda autodock-vina",
        "obabel": "conda install -c conda-forge openbabel",
        "fasttree": "conda install -c bioconda fasttree",
    },
}


def get_tool_version(executable, version_flag):
    """Try to run `executable version_flag` and return a short version string."""
    try:
        cmd = [executable] if not version_flag else [executable, version_flag]
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=10,
        )
        output = (proc.stdout or "") + (proc.stderr or "")
        first_line = output.strip().splitlines()[0] if output.strip() else "(no version output)"
        return first_line[:80]
    except Exception as exc:  # noqa: BLE001
        return f"(found but version check failed: {exc})"


def check_cli_tools():
    system = platform.system()
    install_map = INSTALL_COMMANDS.get(system, INSTALL_COMMANDS["Linux"])
    results = []
    for label, executable, version_flag, fallback in CLI_TOOLS:
        path = shutil.which(executable)
        if path:
            version = get_tool_version(executable, version_flag)
            results.append((label, "PASS", f"{path} -> {version}"))
            continue

        if fallback:
            fallback_exe, fallback_flag = fallback
            fallback_path = shutil.which(fallback_exe)
            if fallback_path:
                version = get_tool_version(fallback_exe, fallback_flag)
                results.append((
                    label,
                    "PASS (substitute)",
                    f"'{executable}' not found, but using {fallback_exe} instead: "
                    f"{fallback_path} -> {version}",
                ))
                continue

        install_cmd = install_map.get(executable, "(no install command known for this OS)")
        results.append((label, "NOT FOUND", f"install with: {install_cmd}"))
    return results, system


# ---------------------------------------------------------------------------
# 3. Report
# ---------------------------------------------------------------------------

def print_table(title, rows):
    print(f"\n{title}")
    print("-" * len(title))
    name_width = max(len(r[0]) for r in rows) + 2
    status_width = max(len(r[1]) for r in rows) + 2
    for name, status, detail in rows:
        print(f"{name:<{name_width}}{status:<{status_width}}{detail}")


def main():
    print("=" * 70)
    print("DHFR Pipeline — Environment Check")
    print("=" * 70)
    print(f"Python: {sys.version.split()[0]}  ({sys.executable})")
    print(f"OS: {platform.system()} {platform.release()} ({platform.machine()})")

    package_results = check_packages()
    print_table("Python package imports", package_results)

    cli_results, system = check_cli_tools()
    print_table("Command-line tools", cli_results)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    all_results = package_results + cli_results
    header = f"{'Check':<20}{'Status':<12}"
    print(header)
    print("-" * len(header))
    n_pass = 0
    n_fail = 0
    for name, status, _ in all_results:
        print(f"{name:<20}{status:<12}")
        if status.startswith("PASS"):
            n_pass += 1
        else:
            n_fail += 1

    print("-" * len(header))
    print(f"{n_pass} PASS / {n_fail} NOT FOUND or FAIL out of {len(all_results)} checks")

    missing_tools = [name for name, status, _ in cli_results if status == "NOT FOUND"]
    if missing_tools:
        print(
            f"\nNote: missing CLI tools ({', '.join(missing_tools)}) must be installed "
            "on the machine actually running the pipeline (this script cannot install "
            "them itself). See the install commands in the table above for this OS "
            f"({system})."
        )

    failed_packages = [name for name, status, _ in package_results if status == "FAIL"]
    if failed_packages:
        print(
            f"\nNote: package import failures ({', '.join(failed_packages)}) usually mean "
            "requirements.txt was not installed into the Python environment currently "
            "running this script. Re-run: pip install -r requirements.txt"
        )

    sys.exit(1 if n_fail > 0 else 0)


if __name__ == "__main__":
    main()
