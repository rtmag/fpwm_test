# 1) Open Project
devtools::document()
devtools::build()

# 2) create and init proj
git init
git remote add origin https://github.com/rtmag/FPWM.git
git add .
git commit -m "FPWM initial commit."
git pull
git push
