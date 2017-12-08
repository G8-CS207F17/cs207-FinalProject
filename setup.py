from setuptools import setup

setup(name='chemkin8',
      version='1.2',
      description='A Chemical Kinetics Library',
      url='https://github.com/G8-CS207F17/cs207-FinalProject',
      author='Mehul Smriti Raje, Ruiqi Chen, Ziqi Guo,',
      author_email='zguo@g.harvard.edu',
      license='Harvard',
      packages=['chemkin8'],
      install_requires=['numpy'],
      setup_requires=['pytest-runner'],
      test_suite='pytest',
      tests_require=['pytest'],
      package_data={'chemkin8':['*', 'tests/*']},
      zip_safe=False)
