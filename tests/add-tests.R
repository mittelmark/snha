library(Rpkg)

x = 2 + 2

if (x != 4) {
    stop("Error: Something strange happened!!")
}

x = add(2,2)

if (x != 4) {
    stop("Error: Something strange happened in add(2,2)!!")
}
