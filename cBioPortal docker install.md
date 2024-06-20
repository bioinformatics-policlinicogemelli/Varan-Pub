
### cBioPortal docker-compose

⚠️ <ins>If cBioPortal is already locally installed, you can skip this paragraph.</ins><br>
⚠️ The operations stated below refer to cBioportal version 5.3.8

This section of the guide will explain how to install the cBioPortal docker image  ([Docker Compose](https://github.com/cBioPortal/cbioportal-docker-compose)) and will show how to load a study folder created with varan 2.0 and how to delete it

To correctly install and use the cBioPortal docker image you can check the following steps:

1. Download [Docker](https://www.docker.com/) to the official website.                                                                                                                      ⚠️ if you use a mac device (Apple Silicon m1,m2) go to the docker settings, features in development and activate the option "Use Rosetta for x86/amd64 emulation on Apple Silicon"

2. Open a terminal and digit the commands:

```
git clone https://github.com/cBioPortal/cbioportal-docker-compose.git
cd cbioportal-docker-compose
./init.sh
```

3. Open the hidden file .env (💡 *Tip*: Run `ls -a` command to look for the .env file inside the repository folder)

```
nano .env
```

4. Set the values of the <b>DOCKER_IMAGE_CBIOPORTAL</b> and <b>DOCKER_IMAGE_SESSION_SERVICE</b> lines as the ones below: \

`DOCKER_IMAGE_CBIOPORTAL=cbioportal/cbioportal:5.3.8` \

`DOCKER_IMAGE_SESSION_SERVICE=cbioportal/session-service:0.6.1`

5. Check that the cBioPortal version of the .env file also matches the ones in <b>docker-compose.arm64.yml </b> and <b>docker-compose.override.yml</b> files (if present).

6. Launch the command:

```
docker-compose up
```
Here are some examples of how to import and delete studies from the local cBioPortal instance:

#### *Import Study*

```
sudo docker compose run cbioportal metaImport.py -u http://cbioportal:8080 -s study/name_example_study/ -o
sudo docker compose restart cbioportal
```

#### *Delete Study*

```
docker compose run cbioportal cbioportalImporter.py -c remove-study -id name_example_study
```
⚠️ For more info about cBioPortal Docker, check [here](https://docs.cbioportal.org/deployment/docker/)
