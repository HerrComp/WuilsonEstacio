// para enlazar con los pc de la unal fisica
ssh -p 443 username@168.176.35.111 
//luego debes ingresar tu contraseña dos veces, el nombre de usauario es la parte principal de correo. 
username: waestacioro
// esto se hace para poder exportar 
export https_proxy="http://username:password@proxyapp.unal.edu.co:8080/"

// luego de esto ya se puede hacer git pull y git push

git push origin master

//para saber con cuantos procesadores estas traajando se hace 

cat /proc/cpuinfo | grep process | wc -l
