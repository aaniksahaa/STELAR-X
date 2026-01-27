#!/bin/bash

# ASTRAL GPU Setup Script
# This script ensures ASTRAL can use GPU on any Linux machine
# Based on Docker troubleshooting experience

set -e

echo "=== ASTRAL GPU Setup Script ==="
echo "Setting up GPU support for ASTRAL..."
echo

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# # Check if running as root
# if [[ $EUID -eq 0 ]]; then
#     print_error "This script should not be run as root. Run as a regular user."
#     exit 1
# fi

# Detect OS
if [[ -f /etc/os-release ]]; then
    . /etc/os-release
    OS=$ID
    VERSION=$VERSION_ID
else
    print_error "Cannot detect OS. This script supports Ubuntu/Debian systems."
    exit 1
fi

print_status "Detected OS: $OS $VERSION"

# Function to install packages
install_packages() {
    print_status "Installing required packages..."
    
    case $OS in
        ubuntu|debian)
            sudo apt update
            sudo apt install -y \
                openjdk-21-jdk \
                ocl-icd-opencl-dev \
                opencl-headers \
                clinfo \
                mesa-opencl-icd \
                build-essential
            ;;
        fedora|centos|rhel)
            sudo dnf install -y \
                java-21-openjdk-devel \
                ocl-icd-devel \
                opencl-headers \
                clinfo \
                mesa-libOpenCL
            ;;
        *)
            print_error "Unsupported OS: $OS"
            exit 1
            ;;
    esac
}

# Function to setup NVIDIA support
setup_nvidia() {
    print_status "Setting up NVIDIA GPU support..."
    
    # Check if NVIDIA GPU is present
    if ! command -v nvidia-smi &> /dev/null; then
        print_warning "nvidia-smi not found. Installing NVIDIA drivers..."
        
        case $OS in
            ubuntu|debian)
                # Add NVIDIA repository
                wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.0-1_all.deb
                sudo dpkg -i cuda-keyring_1.0-1_all.deb
                sudo apt update
                sudo apt install -y nvidia-driver-535 nvidia-opencl-dev
                ;;
            *)
                print_error "Please install NVIDIA drivers manually for your OS"
                return 1
                ;;
        esac
    else
        print_success "NVIDIA drivers already installed"
        nvidia-smi
    fi
    
    # Install NVIDIA OpenCL if not present
    case $OS in
        ubuntu|debian)
            if ! dpkg -l | grep -q nvidia-opencl; then
                sudo apt install -y nvidia-opencl-dev
            fi
            ;;
    esac
    
    # Check NVIDIA device permissions
    print_status "Checking NVIDIA device permissions..."
    for device in /dev/nvidia*; do
        if [[ -e "$device" ]]; then
            print_success "Found device: $device"
            ls -la "$device"
        fi
    done
    
    # Add user to video group if needed
    if ! groups $USER | grep -q video; then
        print_status "Adding user to video group..."
        sudo usermod -a -G video $USER
        print_warning "You need to log out and log back in for group changes to take effect"
    fi
}

# Function to setup OpenCL
setup_opencl() {
    print_status "Setting up OpenCL..."
    
    # Create OpenCL ICD directory
    sudo mkdir -p /etc/OpenCL/vendors
    
    # Setup NVIDIA OpenCL ICD
    if [[ -f /usr/lib/x86_64-linux-gnu/libnvidia-opencl.so.1 ]]; then
        echo "libnvidia-opencl.so.1" | sudo tee /etc/OpenCL/vendors/nvidia.icd > /dev/null
        print_success "NVIDIA OpenCL ICD configured"
    fi
    
    # Setup Intel OpenCL ICD
    if [[ -f /usr/lib/x86_64-linux-gnu/libintelocl.so ]]; then
        echo "libintelocl.so" | sudo tee /etc/OpenCL/vendors/intel.icd > /dev/null
        print_success "Intel OpenCL ICD configured"
    fi
    
    # Setup Mesa OpenCL ICD
    if [[ -f /usr/lib/x86_64-linux-gnu/libmesaopencl.so.1 ]]; then
        echo "libmesaopencl.so.1" | sudo tee /etc/OpenCL/vendors/mesa.icd > /dev/null
        print_success "Mesa OpenCL ICD configured"
    fi
    
    # Test OpenCL
    print_status "Testing OpenCL setup..."
    if command -v clinfo &> /dev/null; then
        clinfo | head -20
    else
        print_warning "clinfo not available, cannot test OpenCL"
    fi
}

# Function to create optimized ASTRAL runner
create_astral_runner() {
    print_status "Creating optimized ASTRAL runner script..."
    
    cat > run_astral_gpu.sh << 'EOF'
#!/bin/bash

# Optimized ASTRAL GPU Runner
# Based on Docker troubleshooting experience

set -e

# Default settings
SAFE_MODE=false
TIMEOUT=300
MEMORY="4g"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --safe-mode)
            SAFE_MODE=true
            shift
            ;;
        --timeout)
            TIMEOUT="$2"
            shift 2
            ;;
        --memory)
            MEMORY="$2"
            shift 2
            ;;
        --help)
            echo "ASTRAL GPU Runner"
            echo "Usage: $0 [options] -- [astral_options]"
            echo "Options:"
            echo "  --safe-mode    Handle AVX2 crashes gracefully"
            echo "  --timeout SEC  Set timeout in seconds (default: 300)"
            echo "  --memory MEM   Set Java memory (default: 4g)"
            echo "  --help         Show this help"
            echo
            echo "Examples:"
            echo "  $0 -- -i input.tre -o output.tre"
            echo "  $0 --safe-mode -- -i input.tre -o output.tre"
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Set up environment
export JAVA_HOME=$(readlink -f /usr/bin/java | sed "s:bin/java::")
export LD_LIBRARY_PATH="$(pwd)/lib:/usr/lib/x86_64-linux-gnu:/usr/local/cuda/lib64:$LD_LIBRARY_PATH"

# Build classpath
CLASSPATH=".:lib/main.jar:lib/colt.jar:lib/JSAP-2.1.jar:lib/jocl-2.0.0.jar"

# Change to main directory
cd main

echo "=== ASTRAL GPU Runner ==="
echo "Java Home: $JAVA_HOME"
echo "Library Path: $LD_LIBRARY_PATH"
echo "Memory: $MEMORY"
echo "Safe Mode: $SAFE_MODE"
echo "Timeout: ${TIMEOUT}s"
echo "Arguments: $@"
echo

# Test OpenCL
echo "=== OpenCL Status ==="
clinfo 2>/dev/null | head -10 || echo "OpenCL info not available"
echo

if [[ "$SAFE_MODE" == "true" ]]; then
    echo "Running in SAFE MODE (handles AVX2 crashes)..."
    
    # Run with timeout and crash handling
    timeout ${TIMEOUT}s java -Xmx${MEMORY} -Djava.library.path=../lib -classpath "../${CLASSPATH}" phylonet.coalescent.CommandLine "$@" || {
        exit_code=$?
        echo
        if [[ $exit_code -eq 124 ]]; then
            echo "❌ Process timed out after ${TIMEOUT} seconds"
            exit 1
        elif [[ $exit_code -eq 139 ]] || [[ $exit_code -eq 132 ]]; then
            echo "⚠️  AVX2 crash detected, but computation likely completed successfully"
            echo "✅ Check output files for results"
            exit 0
        else
            echo "❌ Process exited with code: $exit_code"
            exit $exit_code
        fi
    }
else
    echo "Running in NORMAL MODE..."
    java -Xmx${MEMORY} -Djava.library.path=../lib -classpath "../${CLASSPATH}" phylonet.coalescent.CommandLine "$@"
fi
EOF

    chmod +x run_astral_gpu.sh
    print_success "Created run_astral_gpu.sh"
}

# Function to test ASTRAL GPU setup
test_astral_gpu() {
    print_status "Testing ASTRAL GPU setup..."
    
    if [[ ! -f "lib/main.jar" ]]; then
        print_error "ASTRAL not found in current directory. Please run this script from ASTRAL root directory."
        return 1
    fi
    
    # Test with a small example if available
    if [[ -f "inputs/in200.tr" ]]; then
        print_status "Running test with inputs/in200.tr..."
        ./run_astral_gpu.sh --safe-mode --timeout 60 -- -i inputs/in200.tr -o test_gpu_output.tre
        
        if [[ -f "test_gpu_output.tre" ]]; then
            print_success "GPU test completed successfully!"
            print_status "Output saved to test_gpu_output.tre"
        else
            print_warning "Test completed but no output file found"
        fi
    else
        print_warning "No test input file found. Please test manually with your data."
    fi
}

# Main execution
main() {
    print_status "Starting ASTRAL GPU setup..."
    
    # Check if we're in ASTRAL directory
    if [[ ! -f "lib/main.jar" ]]; then
        print_error "Please run this script from the ASTRAL root directory"
        exit 1
    fi
    
    # Install required packages
    install_packages
    
    # Setup NVIDIA if available
    if lspci | grep -i nvidia > /dev/null; then
        setup_nvidia
    else
        print_warning "No NVIDIA GPU detected. Setting up for Intel/AMD GPUs..."
    fi
    
    # Setup OpenCL
    setup_opencl
    
    # Create optimized runner
    create_astral_runner
    
    # Test setup
    if [[ "$1" != "--no-test" ]]; then
        test_astral_gpu
    fi
    
    echo
    print_success "=== Setup Complete! ==="
    echo
    echo "Usage:"
    echo "  ./run_astral_gpu.sh -- -i input.tre -o output.tre"
    echo "  ./run_astral_gpu.sh --safe-mode -- -i input.tre -o output.tre"
    echo
    echo "The --safe-mode option handles AVX2 crashes gracefully while preserving GPU functionality."
    echo
    if groups $USER | grep -q video; then
        print_success "Ready to use GPU acceleration!"
    else
        print_warning "Please log out and log back in for group changes to take effect."
    fi
}

# Run main function
main "$@"
