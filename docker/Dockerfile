FROM julia:1.11-bookworm

ENV JULIA_DEPOT_PATH=/opt/julia-depot

ENV DEBIAN_FRONTEND=noninteractive

RUN rm /var/lib/dpkg/info/libc-bin.* \
    && apt-get clean \
    && apt-get update \
    && apt-get install -y \
        git \
        libcurl4-openssl-dev \
        cmake \
        wget \
        procps \
        clang \
        libblas-dev \
        liblapack-dev \
        libopenblas-dev \
        tabix \
        unzip \
        openjdk-17-jdk \
        zlib1g-dev

# Install PLINK
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20250819.zip \
    && unzip plink_linux_x86_64_20250819.zip \
    && mv plink /usr/local/bin \
    && chmod +x /usr/local/bin/plink
    
# Install PLINK2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20250129.zip \
    && unzip plink2_linux_x86_64_20250129.zip \
    && mv plink2 /usr/local/bin \
    && chmod +x /usr/local/bin/plink2

COPY assets/ /opt/PgenInteractions/assets/

COPY bin/ /opt/PgenInteractions/bin/

COPY Project.toml /opt/PgenInteractions/Project.toml

COPY src/ /opt/PgenInteractions/src/

RUN julia --project=/opt/PgenInteractions --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' \
    && julia --project=/opt/PgenInteractions --startup-file=no --threads=auto -e 'using PackageCompiler; create_sysimage(["PgenInteractions"], sysimage_path="/opt/PgenInteractions/sysimage.so", include_transitive_dependencies=false, cpu_target=PackageCompiler.default_app_cpu_target())'