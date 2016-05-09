from setuptools import setup


import versioneer

versioneer.versionfile_source = 'multilevelfiredrakeproject/_version.py'
versioneer.versionfile_build = 'multilevelfiredrakeproject/_version.py'
versioneer.tag_prefix = 'v'
versioneer.parentdir_prefix = 'multilevelfiredrakeproject-'
versioneer.VCS = "git"



cmdclass = versioneer.get_cmdclass()

setup(name='multilevelfiredrake',
      version=versioneer.get_version(),
      description='multilevelfiredrake is a package developed to allow the construction of multilevel monte carlo estimators for the uncertainty quantification in firedrake-solvable PDE problems with random components',
      url='http://github.com/alsgregory/multilevelfiredrakeproject',
      author='Alastair C. A. Gregory',
      author_email='alsgregory@aol.com',
      license='Imperial College London',
      packages=['multilevelfiredrakeproject'],
      zip_safe=False)
