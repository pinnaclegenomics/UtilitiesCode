#!/bin/bash

INSTALLDIR=~/bin

echo " "
echo "- Installing utility.scripts..."
echo " "

echo "1. Creating "${INSTALLDIR}" if it does not exist..."
mkdir -p ${INSTALLDIR}

echo "2. Adding "${INSTALLDIR}" to last line of .bashrc/.profile, if it is not already there..."
case :$PATH: in
  *:${INSTALLDIR}:*)
	  echo "   "${INSTALLDIR}" already in PATH"
	  ;;  
  *)
	  FILE1=~/.bashrc
	  FILE2=~/.profile
	  if ! [ -f "$FILE1" -o -f "$FILE2" ]; then
		  #both .bashrc and .profile are missing, so just create a .bashrc
		  echo "export PATH=${INSTALLDIR}:"'$PATH' >> ~/.bashrc
		  echo "   export PATH=${INSTALLDIR}:"'$PATH'" line added to .bashrc"
		  echo "   type source ~/.profile to reload the .bashrc file"
	  elif [ -f "$FILE1" ]; then
		  #check if the line to be added is already in the .bashrc, and if not, add it
		  if $(grep "export PATH=${INSTALLDIR}:"'$PATH' ~/.bashrc); then
			  echo "   export PATH=${INSTALLDIR}:"'$PATH'" already in .bashrc"
		  else
			  echo "export PATH=${INSTALLDIR}:"'$PATH' >> ~/.bashrc
			  echo "   export PATH=${INSTALLDIR}:"'$PATH'" line added to .bashrc"
			  echo "   type source ~/.profile to reload the .bashrc file"
		  fi   
	  elif [ -f "$FILE2" ]; then
		  #check if the line to be added is already in the .profile, and if not, add it
		  if $(grep "export PATH=${INSTALLDIR}:"'$PATH' ~/.profile); then
			  echo "   export PATH=${INSTALLDIR}:"'$PATH'" already in .profile"
		  else
			  echo "export PATH=${INSTALLDIR}:"'$PATH' >> ~/.profile
			  echo "   export PATH=${INSTALLDIR}:"'$PATH'" line added to .profile"
			  echo "   type source ~/.profile to reload the .profile file"
		  fi
	  fi
	  ;;
esac

echo "3. Copying scripts to "${INSTALLDIR}"..."
cp mantaVcfToBedpe/mantaVcfToBedpe.R ${INSTALLDIR}/mantaVcfToBedpe
cp strelkaVcfFilter/strelkaVcfFilter.R ${INSTALLDIR}/strelkaVcfFilter
cp cavemanVcfFilter/cavemanVcfFilter.R ${INSTALLDIR}/cavemanVcfFilter
cp pindelVcfFilter/pindelVcfFilter.R ${INSTALLDIR}/pindelVcfFilter
cp brassBedpeFilter/brassBedpeFilter.R ${INSTALLDIR}/brassBedpeFilter
cp ascatToTsv/ascatToTsv.R ${INSTALLDIR}/ascatToTsv
chmod -R 755 ${INSTALLDIR}
echo "   done."

echo " "
