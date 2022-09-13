# -*- coding: utf-8 -*-
from setuptools import setup

package_dir = \
{'': 'src'}

packages = \
['hapsmash']

package_data = \
{'': ['*']}

install_requires = \
['argparse>=1.4.0,<2.0.0',
 'natsort>=8.0.2,<9.0.0',
 'scipy>=1.7.3,<2.0.0',
 'svglib>=1.4.1,<2.0.0',
 'svgwrite>=1.4.3,<2.0.0']

extras_require = \
{':python_version >= "3.8" and python_full_version < "4.0.0"': ['numpy>=1.22.1,<2.0.0']}

entry_points = \
{'console_scripts': ['hapsmash = hapsmash.__main__:main']}

setup_kwargs = {
    'name': 'hapsmash',
    'version': '0.0.1',
    'description': 'Hapsmash',
    'long_description': '',
    'author': 'Sangjin Lee',
    'author_email': 'sl17@sanger.ac.uk',
    'maintainer': None,
    'maintainer_email': None,
    'url': 'https://github.com/sjin09/hapsmash',
    'package_dir': package_dir,
    'packages': packages,
    'package_data': package_data,
    'install_requires': install_requires,
    'extras_require': extras_require,
    'entry_points': entry_points,
    'python_requires': '>=3.7,<3.11',
}


setup(**setup_kwargs)

