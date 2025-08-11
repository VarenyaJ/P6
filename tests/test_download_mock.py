"""
Isolated test for the 'download' command without hitting the network.

We patch requests.get twice:
- first for the 'latest release' lookup (returns {'tag_name': 'vX'})
- second for the actual file download (returns the file content).
"""

from click.testing import CliRunner
from unittest.mock import patch, Mock
from P6.__main__ import main


def test_download_mocks_network(tmp_path):
    runner = CliRunner()

    def fake_get(url, *args, **kwargs):
        if url.endswith("/releases/latest"):
            return Mock(status_code=200, json=lambda: {"tag_name": "vX"})
        # second call returns the content of hp.json
        return Mock(status_code=200, content=b"{}")

    with patch("P6.__main__.requests.get", side_effect=fake_get):
        res = runner.invoke(main, ["download", "-d", str(tmp_path)])
        assert res.exit_code == 0
        assert (tmp_path / "hp.json").exists()
