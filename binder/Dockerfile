FROM kevincoakley/metenp-base:20220816

ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

# Set the working directory to the home directory
WORKDIR ${HOME}

# Copy the repository to the home directory
COPY . ${HOME}

USER root

# Update the permissions of files in the home directory
RUN chown -R ${NB_UID}:${NB_UID} ${HOME}

USER ${USER}

# Install metenp
RUN R -e 'devtools::install(".", upgrade = "never")'