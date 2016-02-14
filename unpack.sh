#!/bin/bash
DATE=$(date)
echo $DATE > unpack.log
mkdir unpacked_data
for i in *.seed
	do
	echo Unpacking $i >> unpack.log
	rdseed -d -o 1 -f $i

	# merges files
	year=2000
	while [ $year -le 2016 ]
		do
		days=0
		while [ $days -le 4 ]
			do
			FileNum=0
			k=0
			mkdir temp
			
			# copies over batches of 100 days to /temp/
			for j in $year.$days*.SAC
				do
				if [ -e $j ]
					then
					mv $j temp/
					echo Copying $j >> unpack.log
					if [ $FileNum = 0 ]
						then
						k=$j
					fi
					(( FileNum++ ))
				fi
			done
			
			X1files=`ls temp | grep X1`
			for m in $X1files
				do
				if [ -e temp/$m ]
					then
					mkdir temp/tempX1
					mv temp/$m temp/tempX1
					sac decimateX1.sm
					echo $m decimated by X1 >> unpack.log
					mv temp/tempX1/$m temp/
					rm -r temp/tempX1
				fi
			done

			XJfiles=`ls temp | grep XJ`
			for m in $XJfiles
				do
				if [ -e temp/$m ]
					then
					mkdir temp/tempXJ
					mv temp/$m temp/tempXJ
					sac decimateXJ.sm
					echo $m decimated by XJ >> unpack.log
					mv temp/tempXJ/$m temp/
					rm -r temp/tempXJ
				fi
			done



			# if more than one file in /temp/, merges them
			if [ $k != 0 ]
				then
				if [ $FileNum != 1 ]
					then
		       			sac merge.sm
				fi
				echo $k is the first >> unpack.log
				mv temp/$k unpacked_data/
				k=0
			fi
			rm -r temp/
			days=`echo $[$days+1]`
		done
		year=`echo $[$year+1]`
	done	
done
