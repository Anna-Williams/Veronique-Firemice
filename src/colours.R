library(pals) # for palettes with large n #kelly()22, #polychrome()#36, cols

# remove the black and white from the pallete, and the light blue similars to cols25
# assessed with pal.bands and pal.cube
kelly_col <- unname(kelly()[-c(1,2,6)])
# remove the orange colour that is similar to the kelly, and sort with prettiest colours first and ugliest or still similar to kelly last
cols24 <- unname(cols25()[c(19,22:24,8:14, 1:4, 15:17,6,7,25,18,20,21)])
# merge all pallettes for long list colours, 
# my favourite are the first 22 (kelly), cols24 is not bad and I tried to keep it distinct
# and polychrome is just added at the end in case we are missing levels (no checked it is safe with the other 42)
cols <- c(kelly_col, cols24, unname(polychrome()))
#Cols bellow checked with pal.safe
col_wt_ko <- c("#666666", "#E25822")

col_magenta_green <- c("#D83472", "#008000") # #D83472 ##df5697 #6eb76e
#test_col_magenta_green <- c("#6eb76e", "#008000" , "#df5697", "#ce0665", "#D83472", "#6eb76e")
