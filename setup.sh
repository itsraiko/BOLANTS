#!/bin/bash
# BOLANTS setup script
set -e

echo "Installing Boltz-2 model files..."
python setup_boltzina.py

echo "Installing MAXIT..."
if [ ! -f "maxit-v11.300-prod-src.tar.gz" ]; then
    wget https://sw-tools.rcsb.org/apps/MAXIT/maxit-v11.300-prod-src.tar.gz
fi
if [ ! -d "maxit-v11.300-prod-src" ]; then
    tar xvf maxit-v11.300-prod-src.tar.gz
fi
if [ ! -f "maxit-v11.300-prod-src/bin/maxit" ]; then
    cd maxit-v11.300-prod-src && make binary && cd ..
fi

echo ""
echo "========================================================"
echo "Setup complete!"
echo ""
echo "Add the following to your ~/.bashrc or ~/.zshrc:"
echo "  export RCSBROOT=\$(realpath maxit-v11.300-prod-src)"
echo "  export PATH=\$PATH:\$RCSBROOT/bin"
echo ""
echo "PLANTS must be installed separately:"
echo "  https://github.com/purnawanpp/plants"
echo "  chmod +x /path/to/PLANTS"
echo "  Set plants_bin in your config.json"
