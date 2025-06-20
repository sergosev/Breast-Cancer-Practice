#!/bin/bash

echo "📦 Setting up virtual environment for Breast Cancer DEA project"

# 1. Check if venv exists
if [ ! -d "../env" ]; then
    echo "🔧 Creating virtual environment..."
    python3 -m venv ../env
fi

# 2. Activate it
source ../env/bin/activate
echo "✅ Virtual environment activated."

# 3. Install packages
REQ_FILE="./requirements.txt"
if [ -f "$REQ_FILE" ]; then
    echo "📄 Installing from $REQ_FILE..."
    pip install --upgrade pip
    pip install -r "$REQ_FILE"
else
    echo "❌ requirements.txt not found at $REQ_FILE. Aborting."
    deactivate
    exit 1
fi

echo "🎉 Setup complete! You can now run your scripts inside the environment."

