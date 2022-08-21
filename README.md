### statistical learning approaches for the study of plant resistance to pathogens
Internship in AgroParisTech 2022 

[Internship Report](https://www.overleaf.com/read/njgjqhvygjcp)

#### Part I - Improved Eagle Algorithm

#### Part II - Hybrid Model

#### Part III - Two-Way Mixed Effects Model

#### Part IIII - State-of-the-Art in Deep Learning

###### Managing environment with conda (Bash)
```
$ conda activate base
(base)$ conda create -n pytorch_env
(base)$ conda activate pytorch_env
(pytorch_env)$ conda install jupyterlab cudatoolkit pytorch scikit-learn
(pytorch_env)$ jupyter-lab --generate-config
```
###### JupyterLab Configuration (Server: ~/.jupyter/jupyter_lab_config.py)
```
(pytorch_env)$ emacs -nw ~/.jupyter/jupyter_lab_config.py

c.ServerApp.ip = '0.0.0.0'
c.ServerApp.notebook_dir = u'/home/myaccount/work/JupyterHome'
c.ServerApp.open_browser = False
c.ServerApp.password = u'password'
c.ServerApp.password_required = True
c.ServerApp.port = 8889

(pytorch_env)$ jupyter-lab 
or 
(pytorch_env)$ jupyter-lab --config .jupyter/jupyter_lab_config.py
```

###### SSH Port Forwarding (Client)
