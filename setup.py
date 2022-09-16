from setuptools import setup

setup(
    name='UniProtClient',
    version='0.0.1.0',
    author='Christian W. Feldmann',
    license="MIT",
    packages=['UniProtClient'],
    author_email='christian.w.feldmann@gmail.com',
    description='Convenient access to UniProt REST API.',
    install_requires=['pandas', 'numpy', 'tqdm']
)