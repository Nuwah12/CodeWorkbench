########## Script for running the comparitive compartment code from HOMER
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <path to tag directory 1> <path to tag directory 2> -r/--resolution [resolution] -s/--superResolution [super resolution] -c/--numcores [cpu]"
    exit 1
fi

dir1=$1
dir2=$2
res=50000
sr=100000
c=16

while [[ $# -gt 0 ]]; do
    case "$1" in
	$dir1)
	    shift 2
	    ;;
	$dir2) 
	    shift 2
	    ;;
        -r|--resolution)
            res="$2"
            shift 2
            ;;
        -s|--superResolution)
            sr="$2"
            shift 2
            ;;
	-c|--numcores)
            c="$2"
	    shift 2
	    ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Run homer getHiCcorrDiff.pl
getHiCcorrDiff.pl ${dir1} ${dir2} -res ${res} -superRes ${sr} -cpu ${c}
