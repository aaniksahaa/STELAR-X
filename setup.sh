#!/bin/bash

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}Setting up STELAR-MP environment...${NC}"

# Check if Java 17 is installed
if ! command -v java &> /dev/null; then
    echo -e "${YELLOW}Installing Java 17...${NC}"
    sudo apt update
    sudo apt install -y openjdk-17-jdk
fi

# Set JAVA_HOME
JAVA_HOME=$(dirname $(dirname $(readlink -f $(which java))))
echo "export JAVA_HOME=$JAVA_HOME" >> ~/.bashrc
echo "export PATH=\$JAVA_HOME/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc

# Verify Java version
echo -e "${YELLOW}Java version:${NC}"
java -version

# Check if Maven is installed
if ! command -v mvn &> /dev/null; then
    echo -e "${YELLOW}Installing Maven...${NC}"
    sudo apt install -y maven
fi

# Verify Maven version
echo -e "${YELLOW}Maven version:${NC}"
mvn -version

# Make scripts executable
chmod +x build.sh run.sh

echo -e "${GREEN}Setup completed!${NC}"
echo -e "You can now run:"
echo -e "1. ${YELLOW}./build.sh${NC} to build the project"
echo -e "2. ${YELLOW}./run.sh${NC} to run with default settings"
echo -e "3. ${YELLOW}./run.sh input.tre output.tre CPU_PARALLEL${NC} to run with specific files" 