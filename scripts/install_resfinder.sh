### Install ResFinder ###

# Install
git clone -b 4.2.0 https://bitbucket.org/genomicepidemiology/resfinder.git
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git resfinder/db_pointfinder
git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git resfinder/db_resfinder
git clone https://bitbucket.org/genomicepidemiology/kma.git resfinder/cge/kma
cd resfinder/cge/kma && make
cd ../..

# Dependencies:
pip3 install biopython
pip3 install tabulate
pip3 install cgecore
