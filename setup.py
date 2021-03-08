from setuptools import setup, find_packages
setup(
  name = 'pykinetic2',
  packages = ['pykinetic2'],
  package_data = {'templates': ['templates/*']},
  install_requires=["numpy", ],
  python_requires='>=3.7',
  include_package_data=True,
)
