for file in `find src`; do 
  if [ ! -f $file ]; then continue; fi; 
  newf="../MCFM-6.8-addon-clean/$file"; 
  echo "*** $file ***"; diff $file $newf; 
  if [ `diff $file $newf | wc -l` -ge 1 ]; then
    echo "cp $newf $file"
  fi
done
