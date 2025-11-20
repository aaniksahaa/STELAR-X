#!/bin/bash

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=== STELAR-MP Build Script ==="

# # Build CUDA code
# echo -e "\n${YELLOW}Building CUDA code...${NC}"
# cd cuda
# make clean
# make
# if [ $? -ne 0 ]; then
#     echo -e "${RED}CUDA compilation failed!${NC}"
#     exit 1
# fi
# cd ..
# echo -e "${GREEN}CUDA compilation successful${NC}"

# Build Java code
echo -e "\n${YELLOW}Building Java code...${NC}"
mvn clean package
if [ $? -ne 0 ]; then
    echo -e "${RED}Java compilation failed!${NC}"
    exit 1
fi
echo -e "${GREEN}Java compilation successful${NC}"

# Create bin directory for compiled classes
echo -e "\n${YELLOW}Creating binary directory...${NC}"
mkdir -p bin

# Get the path to the Maven repository
M2_REPO="$HOME/.m2/repository"

# Compile Java classes with Maven dependencies in classpath
javac -sourcepath src \
      -d bin \
      -cp "target/stelar-mp-1.0-SNAPSHOT.jar:$M2_REPO/net/java/dev/jna/jna/5.13.0/jna-5.13.0.jar:$M2_REPO/net/java/dev/jna/jna-platform/5.13.0/jna-platform-5.13.0.jar" \
      $(find src -name "*.java")

if [ $? -ne 0 ]; then
    echo -e "${RED}Java class compilation failed!${NC}"
    exit 1
fi
echo -e "${GREEN}Java class compilation successful${NC}"

echo -e "\n${GREEN}All builds completed successfully!${NC}" 