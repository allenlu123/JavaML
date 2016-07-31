import org.jblas.DoubleMatrix;
import java.lang.Math;
import org.jblas.MatrixFunctions;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
public class MachineLearning {
	public void print(DoubleMatrix m){
		for(int i = 0;i < m.getRows();i++){
			int j;
			for(j = 0;j < m.getColumns() - 1;j++){
				System.out.print(Double.toString(m.get(i,j)) + "  ");
			}
			System.out.println(Double.toString(m.get(i,j)));
		}
	}
	private DoubleMatrix findClosestCentroids(DoubleMatrix X, DoubleMatrix centroids){
		DoubleMatrix idx = DoubleMatrix.zeros(X.getRows(),1);
		for(int i=0;i<X.getRows();i++){
			DoubleMatrix row = X.getRow(i);
			DoubleMatrix temp_centroids = MatrixFunctions.pow(MatrixFunctions.pow(centroids.subRowVector(row),2).rowSums(),0.5);	
			int min_index = temp_centroids.columnArgmins()[0];
			idx.put(i, (double) min_index);
		}	
		return idx;
	}
	private DoubleMatrix initKRandomMeans(DoubleMatrix X, int K){
		Integer[] arr = new Integer[X.getRows()];
		for(int i=0;i<arr.length;i++){
			arr[i] = i;
		}
		Collections.shuffle(Arrays.asList(arr));
		DoubleMatrix ret = DoubleMatrix.zeros(K,X.getColumns());
		for(int i=0;i < K;i++){
			ret.putRow(i, X.getRow(arr[i]));
		}
		return ret;
	}
	private DoubleMatrix computeCentroids(DoubleMatrix X, DoubleMatrix idx, int K){
		int n = X.getColumns();
		DoubleMatrix centroids = DoubleMatrix.zeros(K,n);
		for(int i=0;i<K;i++){
			int nums = 0;
			DoubleMatrix members = idx.eq(i);
			nums = (int) members.sum();
			int[] indices = members.findIndices();
			DoubleMatrix cumulative_mean = DoubleMatrix.zeros(1,X.getColumns());
			for(int j=0;j<indices.length;j++){
				cumulative_mean = cumulative_mean.addRowVector(X.getRow(indices[j]));
			}
			cumulative_mean = cumulative_mean.div((double) nums);
			centroids.putRow(i, cumulative_mean);
		}
		return centroids;
	}
	public ArrayList<DoubleMatrix> runKMeans(DoubleMatrix X, int K){
		if(K < 1) return null;
		DoubleMatrix centroids = initKRandomMeans(X,K);
		DoubleMatrix idx = DoubleMatrix.zeros(X.getRows(),1);
		for(int it=0;it<100;it++){
			idx = findClosestCentroids(X,centroids);
			centroids = computeCentroids(X,idx,K);
		}
		ArrayList<DoubleMatrix> ret = new ArrayList<DoubleMatrix>();
		ret.add(centroids);
		ret.add(idx);
		return ret;
	}
	private double sumSquaredErrors(DoubleMatrix X,DoubleMatrix centroids, DoubleMatrix idx){
		double sse = 0;
		DoubleMatrix centroid_matrix = DoubleMatrix.zeros(idx.getRows(),centroids.getColumns());
		for(int i=0;i<idx.getRows();i++){
			centroid_matrix.putRow(i, centroids.getRow((int) idx.get(i,0)));
		}
		sse = MatrixFunctions.pow(X.sub(centroid_matrix), 2).sum();
		return sse;
	}
	private double[] gapStatistic(DoubleMatrix X,int K) throws InterruptedException, ExecutionException{
		DoubleMatrix mins = X.columnMins();
		DoubleMatrix maxes = X.columnMaxs();
		DoubleMatrix Bwkbs = DoubleMatrix.zeros(10,1);
		ArrayList<DoubleMatrix> arr = runKMeans(X,K);
		double Wk = Math.log(sumSquaredErrors(X,arr.get(0),arr.get(1)));
		List<DoubleMatrix> Xbs = new ArrayList<DoubleMatrix>();
		for(int i = 0;i < 10;i++){
			DoubleMatrix Xb = DoubleMatrix.rand(X.getRows(),X.getColumns());
			for(int j = 0;j < X.getColumns();j++){
				Xb.putColumn(j,Xb.getColumn(j).mmul(maxes.get(0,j) - mins.get(0,j)).add(mins.get(0,j)));
			}
			Xbs.add(Xb);
		}
		List<ArrayList<DoubleMatrix>> arrs = runAllKMeans(Xbs,K);
		for(int i = 0;i < 10;i++){
			Bwkbs.put(i,0,Math.log(sumSquaredErrors(Xbs.get(i),arrs.get(i).get(0),arrs.get(i).get(1))));
		}
		double Wkbs = Bwkbs.sum()/10.0;
		double sk = Math.sqrt(MatrixFunctions.pow(Bwkbs.sub(Wkbs),2).sum()/10.0) * Math.sqrt(1.0 + 1.0/10.0);
		return new double[]{Wkbs - Wk,sk};
	}
	private List<ArrayList<DoubleMatrix>> runAllKMeans(List<DoubleMatrix> Xs,int K) 
			throws InterruptedException, ExecutionException {
		int threads = Runtime.getRuntime().availableProcessors();
		final int k = K;
		ExecutorService service = Executors.newFixedThreadPool(threads);
		List<Future<ArrayList<DoubleMatrix>>> futures = new ArrayList<Future<ArrayList<DoubleMatrix>>>();
		for (final DoubleMatrix input : Xs) {
			Callable<ArrayList<DoubleMatrix>> callable = new Callable<ArrayList<DoubleMatrix>>() {
	            public ArrayList<DoubleMatrix> call() throws Exception {
	            	ArrayList<DoubleMatrix> output = runKMeans(input,k);
	                return output;
	            }
			};
			futures.add(service.submit(callable));
		}
		service.shutdown();
		List<ArrayList<DoubleMatrix>> outputs = new ArrayList<ArrayList<DoubleMatrix>>();
	    for (Future<ArrayList<DoubleMatrix>> future : futures) {
	        outputs.add(future.get());
	    }
	    return outputs;
	}
	private List<ArrayList<Double>> runAllGaps(List<Integer> Ks,DoubleMatrix X) 
			throws InterruptedException, ExecutionException {
		int threads = Runtime.getRuntime().availableProcessors();
		final DoubleMatrix x = X.dup();
		ExecutorService service = Executors.newFixedThreadPool(threads);
		List<Future<ArrayList<Double>>> futures = new ArrayList<Future<ArrayList<Double>>>();
		for (final int k : Ks) {
			Callable<ArrayList<Double>> callable = new Callable<ArrayList<Double>>() {
	            public ArrayList<Double> call() throws Exception {
	            	double[] temp = gapStatistic(x,k);
	            	ArrayList<Double> output = new ArrayList<Double>();
	            	output.add(temp[0]);
	            	output.add(temp[1]);
	                return output;
	            }
			};
			futures.add(service.submit(callable));
		}
		service.shutdown();
		List<ArrayList<Double>> outputs = new ArrayList<ArrayList<Double>>();
	    for (Future<ArrayList<Double>> future : futures) {
	        outputs.add(future.get());
	    }
	    return outputs;
	}
	public int chooseCluster(DoubleMatrix X, int maxK) throws InterruptedException, ExecutionException{
		double[] gaps = new double[maxK];
		double[] first = gapStatistic(X,3);
		gaps[0] = first[0];
		int cluster = maxK;
		for(int k = 4;k < maxK;k++){
			double[] gapK1 = gapStatistic(X,k+1);
			gaps[k] = gapK1[0];
			if(Double.compare(gaps[k - 1], gaps[k] - gapK1[1]) >= 0){
				cluster = k;
				break;
			}
		}
		return cluster;
	}
	public HashMap<Restaurant,Integer> getClosestUser(ArrayList<DoubleMatrix> arr, DoubleMatrix X){
		DoubleMatrix centroid_indices = arr.get(1);
		int[] other_indices_and_user = centroid_indices.eq(centroid_indices.get(centroid_indices.getRows() - 1)).findIndices();
		Arrays.sort(other_indices_and_user);
		int[] other_indices = Arrays.copyOfRange(other_indices_and_user, 0, other_indices_and_user.length - 1);
		int[] cols = new int[X.getColumns()];
		for(int i = 0;i < cols.length;i++){
			cols[i] = i;
		}
		DoubleMatrix other_users = X.get(other_indices,cols);
		DoubleMatrix sse = MatrixFunctions.pow(other_users.subRowVector(X.getRow(X.getRows() - 1)), 2.0).rowSums();
		DoubleMatrix closest = X.getRow(other_indices[sse.argmin()]);
		AzureConnector az = new AzureConnector();
		ArrayList<Double> means = az.getMeans("MEANS");
		for(int i = 0;i < closest.getColumns();i++){
			if(Double.compare(closest.get(0,i), means.get(i)) == 0){
				closest.put(0, i, 0.0);
			}
		}
		int[] rest_indices = closest.findIndices();
		ArrayList<Integer> param = new ArrayList<Integer>();
		for(int i : rest_indices){
			param.add(i);
		}
		ArrayList<Restaurant> rests = az.getRestaurantsByNumbers("RESTAURANTS", param);
		ArrayList<String> rest_names = new ArrayList<String>();
		for(Restaurant r : rests){
			rest_names.add(r.getName());
		}
		return az.getRestaurantNumbers("RESTAURANTS", rest_names);
	}
}
