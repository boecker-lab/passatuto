package de.unijena.bioinf.statistics;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Map;
import java.util.TreeMap;

import de.unijena.bioinf.decoy.model.MassBank;

public class Utils {

	
	public static void quicksort(double[] main, int[] index) {
	    quicksort(main, index, 0, index.length - 1);
	}

	// quicksort a[left] to a[right]
	public static void quicksort(double[] a, int[] index, int left, int right) {
	    if (right <= left) return;
	    int i = partition(a, index, left, right);
	    quicksort(a, index, left, i-1);
	    quicksort(a, index, i+1, right);
	}

	// partition a[left] to a[right], assumes left < right
	private static int partition(double[] a, int[] index, 
	int left, int right) {
	    int i = left - 1;
	    int j = right;
	    while (true) {
	        while (less(a[++i], a[right]))      // find item on left to swap
	            ;                               // a[right] acts as sentinel
	        while (less(a[right], a[--j]))      // find item on right to swap
	            if (j == left) break;           // don't go out-of-bounds
	        if (i >= j) break;                  // check if pointers cross
	        exch(a, index, i, j);               // swap two elements into place
	    }
	    exch(a, index, i, right);               // swap with partition element
	    return i;
	}

	// is x < y ?
	private static boolean less(double x, double y) {
		double first=Double.isNaN(x)?Double.NEGATIVE_INFINITY:x;
		double second=Double.isNaN(y)?Double.NEGATIVE_INFINITY:y;
	    return (first > second);
	}

	// exchange a[i] and a[j]
	private static void exch(double[] a, int[] index, int i, int j) {
	    double swap = a[i];
	    a[i] = a[j];
	    a[j] = swap;
	    int b = index[i];
	    index[i] = index[j];
	    index[j] = b;
	}
	

	public static TreeMap<String, MassBank> getAllDBFiles(File[] inputFiles, String ser, boolean forceRecalculation, boolean withPeaks)throws Exception {
		TreeMap<String, MassBank> dbFiles=null;
		if(!forceRecalculation){
			try{
				ObjectInputStream ois = new ObjectInputStream(new FileInputStream(ser));
				dbFiles = (TreeMap<String, MassBank>) ois.readObject();
				ois.close();
			}catch(Exception e){
				System.err.println(e);
			}
		}
		if(dbFiles==null){
			System.out.println("calculating database file...");
			dbFiles=new TreeMap<String,MassBank>();
			for(File f:inputFiles){
				getAllDBFiles(f, dbFiles, withPeaks);
			}
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(ser));
			oos.writeObject(dbFiles);
			oos.flush();
			oos.close();
			System.out.println("...finished");
		}

		return dbFiles;
	}

	public static void getAllDBFiles(File dbFolder, Map<String, MassBank> results, boolean withPeaks) throws Exception{
		if(dbFolder.exists()){
			for(File f:dbFolder.listFiles()){
				if(f.isFile()&&f.getName().endsWith(".txt")){
					results.put(SimilarityMatrix.getIDDatasetChargeOrigin(f.getName().replaceAll(".txt","")), new MassBank(f, withPeaks));
					System.out.println(f.getAbsolutePath());
				}else if (f.isDirectory()){
					getAllDBFiles(f,results, withPeaks);
				}
			}
		}
	}
}
