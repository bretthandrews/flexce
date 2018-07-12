# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 15:07:67
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 16:07:68

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--runslow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow"):
        # --runslow given in cli: do not skip slow tests
        return
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
            print('\nSkipped tests marked as "slow". Use --runslow option to run them.\n')
