# @Author: Brett Andrews <andrews>
# @Date:   2018-07-12 09:07:44
# @Last modified by:   andrews
# @Last modified time: 2018-07-12 10:07:65

"""
FILE
    general.py

DESCRIPTION
    General purpose functions for file I/O.
"""

from os.path import join


def _make_sim_path(path, sim_id=None, ext=None):
    """Construct path to simulation output file.

    Args:
        path (str): Output path.
        sim_id (str): Simulation ID number. Default is ``None``.
        ext (str): File extension (either '.pck' or '.txt'). Default is
            ``None``.

    Returns:
        str: file path and name.
    """
    if sim_id is not None:
        ext = ext if ext[0] == '.' else '.' + ext

        assert ext in ['.pck', '.txt'], '``ext`` must be either ".pck" or ".txt".'

        prefix = {'.pck': 'sim',
                  '.txt': 'ab'}

        sim_id = str(sim_id)
        path = join(path, f'sim{sim_id}', f'{prefix[ext]}{sim_id}{ext}')

    return path
