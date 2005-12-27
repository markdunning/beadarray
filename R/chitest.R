"chitest" <-
function(observed, expected){

#Applies standard chi-square formula to a set of expected and observed values

chi = sum((observed-expected)^2/expected)

chi

}

