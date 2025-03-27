from click.testing import CliRunner
from scsilhouette.cli import cli

def test_cli_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "Usage:" in result.output

def test_download_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["download", "--help"])
    assert result.exit_code == 0
    assert "--url" in result.output

def test_compute_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["compute", "--help"])
    assert result.exit_code == 0
    assert "--h5ad-path" in result.output
