import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = 'pykinetic2',
    packages = ['pykinetic2'],
    description = """A python library and command line apps to write 
                     microkinetic models for """,
    keywords = ['compchem', 'microkinetics'],
    author = 'Raúl Pérez-Soto',
    classifiers = ["Programming Language :: Python :: 3",],
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url = 'https://github.com/fmaserasgroup-repo/pyssian-utils', # Update
    python_requires='>=3.7',
    install_requires=['setuptools','pathlib','numpy'],
    include_package_data=True,
    package_data = {'templates': ['templates/*']},
    scripts = ['pykinetic2/pyssian-thermo.py',
               'pykinetic2/pyssian-potential.py'],
)
