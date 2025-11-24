for f in .github/sv/configs/*.json; do
    BN=$(basename "$f")
    BN="${BN%.*}"
    echo "Running snakevis on config: $BN"
    python3 .github/sv/snakevis.py -c .github/sv/snakevis.json -s "$f" -u upload
done