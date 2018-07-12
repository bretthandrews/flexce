# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 09:07:42
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 10:07:79

import pytest

import flexce.fileio.general


class TestGeneral(object):

    @pytest.mark.parametrize('sim_id', ['0', 0])
    @pytest.mark.parametrize('ext', ['pck', '.pck'])
    def test_make_sim_path_pck_with_sim_id(self, sim_id, ext):
        actual = flexce.fileio.general._make_sim_path('mypath', sim_id, ext)
        expected = 'mypath/sim0/sim0.pck'
        assert actual == expected

    @pytest.mark.parametrize('sim_id', ['0', 0])
    @pytest.mark.parametrize('ext', ['txt', '.txt'])
    def test_make_sim_path_txt_with_sim_id(self, sim_id, ext):
        actual = flexce.fileio.general._make_sim_path('mypath', sim_id, ext)
        expected = 'mypath/sim0/ab0.txt'
        assert actual == expected

    def test_make_sim_path_no_sim_id(self):
        actual = flexce.fileio.general._make_sim_path('mypath/sim0/sim0.pck')
        expected = 'mypath/sim0/sim0.pck'
        assert actual == expected
