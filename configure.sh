#!/bin/bash -l
##### Date:04 June 2019
##### configure.sh: set the current path for loading R functions
istr="/path/to"
ostr="$PWD"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/AEM_update_X_beta.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/Create_count_matrix.R"
eval "sed -i -e 's#"$istr"#"$ostr"#g' R/buildCRP.R"
echo "Done!"