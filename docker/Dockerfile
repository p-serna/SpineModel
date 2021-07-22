FROM python:3.7

ENV USER=p-serna REPO=SpineModel BRANCH=resubmission

RUN apt-get update && apt-get -y install vim 
RUN apt-get -y --no-install-recommends install texlive texlive-latex-extra texlive-science dvipng
RUN apt-get -y --no-install-recommends install cm-super
# mx=320;my=256;head -c "$((3*mx*my))" /dev/urandom | convert -depth 8 -size "${mx}x${my}" RGB:- Neuron_persp9.png
# mx=320;my=256;head -c "$((3*mx*my))" /dev/urandom | convert -depth 8 -size "${mx}x${my}" RGB:- spine-inh_sketch_bothm.png

RUN pip install numpy matplotlib scipy cython ipython
RUN pip install pandas==1.2.5
RUN pip install jupyter
RUN pip install neuron=8.0
RUN mkdir neuron
RUN git -C neuron/ clone https://github.com/neuronsimulator/iv

# To avoid caching git repo if it gets updated... in principle
# https://stackoverflow.com/questions/36996046/how-to-prevent-dockerfile-caching-git-clone
ADD https://api.github.com/repos/$USER/$REPO/git/refs/heads/$BRANCH version.json
RUN git clone --branch $BRANCH https://github.com/$USER/$REPO.git

COPY script.sh .
RUN /script.sh

WORKDIR /$REPO
# https://u.group/thinking/how-to-put-jupyter-notebooks-in-a-dockerfile/
CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
#CMD ipython
#CMD /bin/bash

#docker build -t spinemodel .
#docker run -p 8888:8888 spinemodel
