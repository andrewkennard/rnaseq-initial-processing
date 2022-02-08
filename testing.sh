
string="This-is_an_example"
prefix=${string%%_*}
echo ${prefix}
echo ${prefix##*-}

