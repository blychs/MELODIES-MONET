from melodies_monet import driver

import warnings
warnings.filterwarnings('ignore')

an = driver.analysis()
an.control = 'wrfchem.yaml'
an.read_control()

an.open_models()
an.models

an.open_obs()
an.obs

# for obs in an.obs:
#     print(an.obs[obs])
#     print(an.obs[obs].obj.info())

an.pair_data()
for key in an.paired:
    print(an.paired[key])

an.plotting()
