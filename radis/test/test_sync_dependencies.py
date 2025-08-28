# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 13:11:41 2025

@author: erwan
"""
import re

import toml
import yaml


def strip_version(dep):
    """Strip version information from a dependency string.

    We remove:
    - version constraints (e.g. '>=1.0.0')
    - environment markers (e.g. '; sys_platform == "win32"')
    - extras (e.g. '[test]')
    - URL schemes (e.g. 'git+https://
    - python version markers (e.g. 'python_version >= "3.6"')
    """
    dep = re.sub(r"[\s;].*", "", dep)  # Remove everything after a space or semicolon
    dep = re.sub(r"\[.*\]", "", dep)  # Remove extras
    dep = re.sub(r"^(git\+|https?://)", "", dep)  # Remove URL schemes
    dep = re.sub(r"python_version.*", "", dep)  # Remove python version markers
    return dep


def get_pip_deps_from_pyproject():
    """Extract dependencies from pyproject.toml."""
    try:
        data = toml.load("../../pyproject.toml")
        deps = data.get("project", {}).get("dependencies", [])
        return set(strip_version(dep) for dep in deps)
    except Exception as e:
        print(f"Error reading pyproject.toml: {e}")
        return set()


def get_deps_from_conda_env():
    """Extract dependencies from environment.yml."""
    try:
        with open("../../environment.yml", "r") as f:
            data = yaml.safe_load(f)

        conda_deps = set()
        pip_deps = set()

        for dep in data.get("dependencies", []):
            if isinstance(dep, str):
                dep_str = strip_version(dep)
                # Ignore "pip" (as a package), "python"
                if dep_str in ["python"]:
                    continue
                if dep_str.startswith("pip"):
                    continue
                # Exception for Pytables; it's called Tables in pip.
                # ... https://www.pytables.org/usersguide/installation.html
                if dep_str.startswith("pytables"):
                    dep_str = dep_str[2:]
                conda_deps.add(dep_str)  # Strip version info
            elif isinstance(dep, dict) and "pip" in dep:
                pip_deps.update(strip_version(pip_dep) for pip_dep in dep["pip"])

        return conda_deps | pip_deps  # Merge conda and pip dependencies
    except Exception as e:
        print(f"Error reading environment.yml: {e}")
        return set()


def check_consistency():
    """Compare dependencies between pyproject.toml and environment.yml."""
    pip_deps_pyproject = get_pip_deps_from_pyproject()
    conda_deps_env = get_deps_from_conda_env()

    missing_in_env = pip_deps_pyproject - conda_deps_env
    missing_in_pyproject = conda_deps_env - pip_deps_pyproject

    print("\n=== Dependency Consistency Check ===")
    if missing_in_env:
        print("⚠️  Missing in environment.yml:", missing_in_env)
    if missing_in_pyproject:
        print("⚠️  Missing in pyproject.toml:", missing_in_pyproject)
    if not missing_in_env and not missing_in_pyproject:
        print("✅ Dependencies are consistent!")
    return missing_in_env, missing_in_pyproject


def test_consistency():
    """Ensures that dependencies are consistent between pyproject.toml and environment.yml."""

    missing_in_env, missing_in_pyproject = check_consistency()
    assert (
        not missing_in_env
    ), f"Missing dependencies in environment.yml: {missing_in_env}"
    assert (
        not missing_in_pyproject
    ), f"Missing dependencies in pyproject.toml: {missing_in_pyproject}"


if __name__ == "__main__":
    test_consistency()
