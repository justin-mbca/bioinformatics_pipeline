"""
Test PuLP API compatibility with Snakemake.

This test validates that the installed PuLP version has the list_solvers() method
that Snakemake requires, preventing AttributeError: module 'pulp' has no attribute 'list_solvers'.
"""
import pulp


def test_pulp_has_list_solvers():
    """Test that pulp has the list_solvers method required by Snakemake."""
    assert hasattr(pulp, 'list_solvers'), (
        "PuLP is missing the 'list_solvers' method. "
        "Ensure pulp<=2.7.0 is installed for Snakemake compatibility."
    )


def test_pulp_list_solvers_callable():
    """Test that pulp.list_solvers() can be called successfully."""
    try:
        solvers = pulp.list_solvers(onlyAvailable=True)
        assert isinstance(solvers, list), "list_solvers should return a list"
    except AttributeError as e:
        raise AssertionError(
            f"Failed to call pulp.list_solvers(): {e}. "
            "Ensure pulp<=2.7.0 is installed for Snakemake compatibility."
        )


def test_pulp_version():
    """Test that pulp version is 2.7.0 or earlier."""
    version = pulp.__version__
    major, minor = map(int, version.split('.')[:2])
    
    assert (major, minor) <= (2, 7), (
        f"PuLP version {version} may not be compatible with Snakemake. "
        "Ensure pulp<=2.7.0 is installed."
    )
