#!/bin/bash

echo "📦 Setting up virtual environment for Breast Cancer DGEA project"

# 1. Check for Python 3
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install it first."
    exit 1
fi

# 2. Check if venv exists, create if not
VENV_DIR="./env"
if [ ! -d "$VENV_DIR" ]; then
    echo "🔧 Creating virtual environment..."
    python3 -m venv "$VENV_DIR" || { echo "❌ Failed to create virtual environment."; exit 1; }
fi

# 3. Activate the virtual environment
echo "✅ Activating virtual environment..."
source "$VENV_DIR/bin/activate" || { echo "❌ Failed to activate virtual environment."; exit 1; }

# 4. Upgrade pip inside the virtual environment
echo "🔄 Upgrading pip..."
pip3 install --upgrade pip || { echo "❌ Failed to upgrade pip."; deactivate; exit 1; }

# 5. Install packages from requirements.txt
REQ_FILE="./requirements.txt"
if [ -f "$REQ_FILE" ]; then
    echo "📄 Installing from $REQ_FILE..."
    pip3 install -r "$REQ_FILE" || { echo "❌ Failed to install dependencies."; deactivate; exit 1; }
else
    echo "❌ requirements.txt not found at $REQ_FILE. Aborting."
    deactivate
    exit 1
fi

echo "🎉 Setup complete! Run 'source ./env/bin/activate' to start using the environment."
echo "🦀 After finishing run 'deactivate' to exit the environment"
echo "🦀 Activate each time when resuming to work"