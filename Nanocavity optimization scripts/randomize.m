function randVal = randomize(nominal,stddev)

randVal = nominal;
if stddev > 0
    randVal = normrnd(nominal,stddev);
end