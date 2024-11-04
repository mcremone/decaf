if voms-proxy-info --exists; then
    echo "VOMS proxy exists."
    timeleft=$(voms-proxy-info --timeleft)
    echo "Time left for the proxy: ${timeleft} seconds." 

else
    echo "VOMS proxy is not valid or has expired."
fi
