prod_interno = function(x, y) {sum(x * y)}

norma = function(x) {sqrt(prod_interno(x, x))}

corVector = function(x, y) {prod_interno(x, y) / (norma(x) * norma(y))}