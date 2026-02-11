# 1. Use the official Julia image as a base
FROM julia:1.10

# 2. Set the working directory inside the container
WORKDIR /app

# 3. Copy your dependency files first (for caching speed)
COPY Project.toml Manifest.toml ./

# 4. Install and precompile packages
# This makes sure the container doesn't hang compiling for 10 mins every time you run it
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

# 5. Copy your actual code
COPY src/ribosome.jl src/ribosome.jl
# (Or use "COPY . ." to copy everything in the folder)

# 6. Command to run when the container starts
CMD ["julia", "--project=.", "src/ribosome.jl"]
