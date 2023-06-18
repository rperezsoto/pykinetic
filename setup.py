from setuptools import setup, find_packages

__version__ = "0.1.0"

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name = 'pykinetic',
    packages = find_packages(),#['pykinetic',],
    description = """A python library and command line apps to write 
                     microkinetic models for """,
    keywords = ['compchem', 'microkinetics'],
    author = 'Raúl Pérez-Soto',
    classifiers = ["Programming Language :: Python :: 3",],
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url = 'https://github.com/rperezsoto/pykinetic', # Update
    python_requires='>=3.6',
    install_requires=['setuptools','pathlib','numpy'],
    extras_requires=['numpy','scipy'],
    include_package_data=True,
    package_data = {'templates': ['pykinetic/templates/*'],
                    'tests'    : ['pykinetic/tests/*.py']},
    scripts = ['pykinetic/scripts/pykinetic-model.py',
               'pykinetic/scripts/pykinetic-scan.py'],
)
