"""Unit tests for the dependency bootstrap in ``importdependencies.py``.

The install path is exercised by monkeypatching ``_is_available`` and
``install`` so the tests never touch the network and pass regardless of which
optional packages happen to be present in the runner's environment.
"""

import importdependencies as dep


def test_dependencies_mapping_is_well_formed():
    assert dep.DEPENDENCIES, "expected at least one managed dependency"
    for pip_name, import_name in dep.DEPENDENCIES.items():
        assert isinstance(pip_name, str) and pip_name
        assert isinstance(import_name, str) and import_name


def test_is_available_true_for_stdlib():
    assert dep._is_available('os') is True


def test_is_available_false_for_missing_module():
    assert dep._is_available('definitely_not_a_real_module_xyz') is False


def test_checkdep_skips_packages_already_present(monkeypatch):
    calls = []
    monkeypatch.setattr(dep, '_is_available', lambda name: True)
    monkeypatch.setattr(dep, 'install', lambda pkg: calls.append(pkg))
    assert dep.checkdep() == []
    assert calls == []


def test_checkdep_installs_each_missing_package(monkeypatch):
    calls = []
    monkeypatch.setattr(dep, '_is_available', lambda name: False)
    monkeypatch.setattr(dep, 'install', lambda pkg: calls.append(pkg))
    installed = dep.checkdep()
    # Only required packages are reported as installed (gate a restart)...
    assert set(installed) == set(dep.DEPENDENCIES.keys())
    # ...but optional/perf packages are also attempted.
    assert set(calls) == set(dep.DEPENDENCIES.keys()) | set(dep.OPTIONAL_DEPENDENCIES.keys())


def test_checkdep_survives_install_failure(monkeypatch):
    def boom(pkg):
        raise RuntimeError('no network')
    monkeypatch.setattr(dep, '_is_available', lambda name: False)
    monkeypatch.setattr(dep, 'install', boom)
    # A failed install must not propagate, and the package is not reported as
    # successfully installed.
    assert dep.checkdep() == []
