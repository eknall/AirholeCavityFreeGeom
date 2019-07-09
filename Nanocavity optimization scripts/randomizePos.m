function [randX,randY] = randomizePos(nomX,nomY,stddev)

randX = nomX;
randY = nomY;
if stddev > 0
    randX = normrnd(nomX,stddev);
    randY = normrnd(nomY,stddev);
end