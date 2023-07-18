// <GUI menuname="simpleNMR" shortcut="Ctrl+1" tooltip="simple NMR export" icon="C:\Users\vsmw51\AppData\Roaming\Mestrelab Research S.L\MestReNova\scripts\simpleNMRLogo.png" />

function simpleNMR(){


  function isObjectEmpty(obj) {
    return Object.keys(obj).length === 0;
  }  
  
function splitPathAndFilename(path) {
  var separatorIndex = path.lastIndexOf('/');
  var directoryPath = path.substring(0, separatorIndex);
  var filename = path.substring(separatorIndex + 1);

  var dotIndex = filename.lastIndexOf('.');
  var name = filename.substring(0, dotIndex);
  var extension = filename.substring(dotIndex + 1);

  return {
    directoryPath: directoryPath,
    filename: filename,
    name: name,
    extension: extension
  };
}  




	function getActiveMolecule(aDocWin, aMolPlugin) {
		var molec = aMolPlugin.activeMolecule();

		if (molec.isValid()) {
			return molec;
		}
		if (aDocWin.itemCount("Molecule") === 1) {
			molec = new Molecule(aDocWin.item(0, "Molecule"));
			return molec;
		}
		return undefined;
	}		
	
  // iterate through pages and print out what is on the page
  doc = Application.mainWindow.activeDocument;
  
  // return molecule from document
  mol = getActiveMolecule(doc, Application.molecule);
  
  
  var spectra = {};
  if( mol === undefined ){
    print("molecule is undefined");  
  }
  else {
    print("molecule is valid");
    var smiles = mol.generateSMILES()
    spectra["smiles"] = smiles;
  }

  print("doc.pageCount(), ", doc.pageCount());


// Example usage:
  var path = doc.name;
  var result = splitPathAndFilename(path);
  print("result.directoryPath, ", result.directoryPath); // Output: /path/to
  print("result.filename, ", result.filename); // Output: somefile.txt

  for (i = 0, pageCount = doc.pageCount(); i < pageCount; i++) {
    page = doc.page(i);
    print(i, " ", page.itemCount() );
    for (j = 0, itemCount = page.itemCount(); j < itemCount; j++) {
      spec = new NMRSpectrum(page.item(j));
      
      //var spectitle = new String(spec.title);
      //print( i, j, spec.isValid(), spec.type, spec.subtype, spec.experimentType, spec.title, spec.originalFormat);      
      if( spec.isValid() ){
      	 var spectitle = spec.title + "_" + i;
        spectra[spectitle] = {};
        spectra[spectitle]["origin"] = spec.originalFormat;
        spectra[spectitle]["type"] = spec.type;
        spectra[spectitle]["subtype"] = spec.subtype;
        spectra[spectitle]["experimenttype"] = spec.experimentType;
        spectra[spectitle]["experiment"] = spec.getParam("Experiment");
        spectra[spectitle]["pulsesequence"] = spec.getParam("Pulse Sequence");
        spectra[spectitle]["intrument"] = spec.getParam("Instrument");
        spectra[spectitle]["probe"] = spec.getParam("Probe");
        spectra[spectitle]["datafilename"] = spec.getParam("Data File Name");
        spectra[spectitle]["comment"] = spec.getParam("Comment");
        
        
        
        // process multiplets
        
        var multiplets = spec.multiplets();
                
        spectra[spectitle]["multiplets"] = {};
        spectra[spectitle]["multiplets"]["count"] = multiplets.count;
        spectra[spectitle]["multiplets"]["normValue"] = multiplets.normValue;
        
        
        // loop through multiplets and add information
        for( m=0; m<multiplets.count; m++){
          multiplet = multiplets.at(m);
          spectra[spectitle]["multiplets"][m] = {};
          spectra[spectitle]["multiplets"][m]["delta1"] = multiplet.delta;
          spectra[spectitle]["multiplets"][m]["nH"] = multiplet.nH;
          spectra[spectitle]["multiplets"][m]["realH"] = multiplet.realH;
          spectra[spectitle]["multiplets"][m]["integralValue"] = multiplet.integralValue();
          spectra[spectitle]["multiplets"][m]["category"] = multiplet.category;
          spectra[spectitle]["multiplets"][m]["type"] = multiplet.type;
          
          var jlist = multiplet.jList();
          spectra[spectitle]["multiplets"][m]["jlistcount"] = jlist.count;
          
          var jvalslist = [];
          
          for( q=0; q<jlist.count; q++){
            jvalslist.push(jlist.at(q));
          }
          spectra[spectitle]["multiplets"][m]["jvals"] = jvalslist;          
        }
        
        // loop over peaks in spectrum and add information
        var peaks = spec.peaks()
        spectra[spectitle]["peaks"] = {};
        spectra[spectitle]["peaks"]["count"] = peaks.count;
        for( p=0; p<peaks.count; p++ ){
          var pk = peaks.at(p);
          if( pk.type === 0 ){
            spectra[spectitle]["peaks"][p] = {};
            spectra[spectitle]["peaks"][p]["delta1"] = pk.delta(1);
            spectra[spectitle]["peaks"][p]["delta2"] = pk.delta(2);
            spectra[spectitle]["peaks"][p]["intensity"] = pk.intensity;
            spectra[spectitle]["peaks"][p]["type"] = pk.type;
          }
        }
        
        // loop over integrals in spectrum and add information
        var integrals = spec.integrals();
        spectra[spectitle]["integrals"] = {};
        spectra[spectitle]["integrals"]["count"] = integrals.count;  
        spectra[spectitle]["integrals"]["normValue"] = integrals.normValue;
        spectra[spectitle]["integrals"]["type"] = integrals.type;    
        
        for( p=0; p<integrals.count; p++){
          var integral = integrals.at(p);
          spectra[spectitle]["integrals"][p] = {};
          spectra[spectitle]["integrals"][p]["integralValue"] = integral.integralValue();
          spectra[spectitle]["integrals"][p]["rangeMin1"] = integral.rangeMin(1);
          spectra[spectitle]["integrals"][p]["rangeMin2"] = integral.rangeMin(2);
          spectra[spectitle]["integrals"][p]["rangeMax1"] = integral.rangeMax(1);
          spectra[spectitle]["integrals"][p]["rangeMax2"] = integral.rangeMax(2);
          
        }
      }        
//      else{
//        print("spec is not valid");
//      }
    }
  } 
  
  print(JSON.stringify(spectra,null,4));  
  print("isObjectEmpty, ", isObjectEmpty(spectra) );
  
  if( isObjectEmpty(spectra) ){
    print("No spectra found");
  }
  else{
    // save spectra JSON string to file in data directory
    var jsonfilename = result.directoryPath + "/" + result.name + "_mresnova.json";
    print("jsonfilename, ", jsonfilename );
    
    // save spectra json string to file
    var json_spectra_string = JSON.stringify(spectra,null,4);
    
 	fout = new File(jsonfilename);
 	
 	if (fout.open(File.WriteOnly)) {
 	  sout = new TextStream(fout);
 	  sout.writeln(json_spectra_string);
	}
   fout.close(); 	

  }
  
                                                                                                    
}              