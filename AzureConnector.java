import java.sql.*;

import org.jblas.DoubleMatrix;

import java.io.File;
import java.io.IOException;
import java.lang.Math;
import org.jblas.MatrixFunctions;
import org.jblas.ranges.IntervalRange;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import com.microsoft.sqlserver.jdbc.*;
public class AzureConnector {
	private final String connectionUrl = "jdbc:sqlserver://allenlu.database.windows.net:1433;"
			+ "database=Allen;user=t-allu@allenlu;password=IndianapolisColts12!;"
			+ "encrypt=true;trustServerCertificate=false;hostNameInCertificate=*.database.windows.net;loginTimeout=30;";
	
	private void swap(ArrayList<Double> arr, ArrayList<Restaurant> rs, int lo, int hi) {
		Collections.swap(arr, lo, hi);
		Collections.swap(rs, lo, hi);
		return;
	}
	private int partition(ArrayList<Double> arr, ArrayList<Restaurant> rs, int lo, int hi) {
		int left = lo + 1;
		int right = hi;
		swap(arr,rs, lo,lo + (hi - lo)/2);
		while(left < right) {
			if(Double.compare(arr.get(left), arr.get(lo)) <= 0) {
				left++;
			}
			else {
				swap(arr,rs,left,right - 1);
				right--;
			}
		}
		swap(arr,rs,lo,right - 1);
		return right - 1;
	}
	public void quickSort(ArrayList<Double> arr, ArrayList<Restaurant> rs, int lo, int hi) {
		if(lo >= hi - 1) return;
		int pvtindex = partition(arr, rs, lo, hi);
		quickSort(arr, rs, lo, pvtindex);
		quickSort(arr, rs, pvtindex + 1, hi);
	}
	public boolean tableExists(String table_name){
		Connection con = null;
		try {
			con = DriverManager.getConnection(connectionUrl);
			DatabaseMetaData meta = con.getMetaData();
			  ResultSet res = meta.getTables(null, null, table_name, 
			     new String[] {"TABLE"});
			  while(res.next()) return true;
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return false;
	}
	public void createTable(String[] cols, String table_name){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "CREATE TABLE " + table_name + " " +
		                   "(id INTEGER not NULL,  client VARCHAR(40), ";
				for(int i = 0;i < cols.length - 1;i++){
					sql += " " + cols[i] + " FLOAT, ";
				}
				sql += " " + cols[cols.length - 1] + " FLOAT)";
				stmt.executeUpdate(sql);
				System.out.println("Created table in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public void createTableUsers(String table_name){
		Connection con = null;
		Statement stmt = null;
		try{
			con = DriverManager.getConnection(connectionUrl);
			stmt = con.createStatement();
			String sql = "CREATE TABLE " + table_name + " (num INTEGER not NULL)";
			stmt.executeUpdate(sql);
			System.out.println("Created table in given database...");
		}catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}finally{
		      //finally block used to close resources
		      try{
		         if(stmt!=null)
		            con.close();
		      }catch(SQLException se){
		      }// do nothing
		      try{
		         if(con!=null)
		            con.close();
		      }catch(SQLException se){
		         se.printStackTrace();
		      }
		}
	}
	public void insertNumUsers(String table_name, int num){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "INSERT INTO " + table_name + " VALUES(" + Integer.toString(num) + ")";
				stmt.executeUpdate(sql);
				System.out.println("Inserted number in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public void updateNumUsers(String table_name, int new_num){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "UPDATE " + table_name + " SET num=" + Integer.toString(new_num);
				stmt.executeUpdate(sql);
				System.out.println("Inserted number in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public int getNumUsers(String table_name){
		Connection con = null;
		Statement stmt = null;
		int num = 0;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT num FROM " + table_name;
				ResultSet rs = stmt.executeQuery(sql);
				while(rs.next()){
					num = rs.getInt("num");
				}
				System.out.println("Selected user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return num;
	}
	public void createTableMeans(String[] cols, String table_name){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "CREATE TABLE " + table_name + " " +
		                   "(" + "col" + cols[0] + " FLOAT, ";
				for(int i = 1;i < cols.length - 1;i++){
					sql += " col" + cols[i] + " FLOAT, ";
				}
				sql += " col" + cols[cols.length - 1] + " FLOAT)";
				stmt.executeUpdate(sql);
				System.out.println("Created means table in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public void insertTableMeans(ArrayList<Double> arr, String table_name){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "INSERT INTO " + table_name + " " +
		                   "(";
				String values = "VALUES(";
				for(int i = 0;i < arr.size() - 1 ;i++){
					sql += "col" + Integer.toString(i) + ", ";
					values += Double.toString(arr.get(i)) + ", ";
				}
				sql += "col" + Integer.toString(arr.size() - 1)+ ") ";
				values += Double.toString(arr.get(arr.size() - 1)) + ")";
				stmt.executeUpdate(sql + values);
				System.out.println("Inserted means in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public void insertIntoTable(HashMap<String,ArrayList<Double>> hm,String table_name){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "INSERT INTO " + table_name + "(id, client, ";
				ArrayList<Double> arr = hm.get(hm.keySet().iterator().next());
				for(int j = 0;j < arr.size();j++){
					sql += "col" + Integer.toString(j) + ", ";
				}
				sql = sql.substring(0,sql.length() - 2) + ")";
				String values = " VALUES(";
				int i = 0;
				for(String key : hm.keySet()){
					values += Integer.toString(i) + ", " + "'" + key + "'" + ", ";
					ArrayList<Double> vals = hm.get(key);
					for(int j = 0;j < vals.size();j++){
						values += Double.toString(vals.get(j)) + ", ";
					}
					values = values.substring(0,values.length() - 2) + "), (";
					if((i + 1) % 100 == 0 || (i + 1) == hm.keySet().size()) {
						System.out.println(Integer.toString(i));
						stmt.executeUpdate(sql + values.substring(0,values.length() - 3));
						values = " VALUES(";
					}
					i++;
				}
				System.out.println("Inserted all...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public void createLoginTable(String table_name){
		Connection con = null;
		Statement stmt = null;
		try{
			con = DriverManager.getConnection(connectionUrl);
			stmt = con.createStatement();
			String sql = "CREATE TABLE " + table_name + " (username VARCHAR(255), " +
							"password VARCHAR(255))";
			stmt.executeUpdate(sql);
			System.out.println("Created table in given database...");
		}catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}finally{
		      //finally block used to close resources
		      try{
		         if(stmt!=null)
		            con.close();
		      }catch(SQLException se){
		      }// do nothing
		      try{
		         if(con!=null)
		            con.close();
		      }catch(SQLException se){
		         se.printStackTrace();
		      }
		}
	}
	public boolean insertUserLogin(String table_name, String user_name, String password){
		Connection con = null;
		Statement stmt = null;
		if(tableExists(table_name)){
			try{
				con = DriverManager.getConnection(connectionUrl);
				String exists = "SELECT * FROM " + table_name + " WHERE username='" + user_name.replace("'", "") + "'";
				stmt = con.createStatement();
				ResultSet rs = stmt.executeQuery(exists);
				while(rs.next()){
					con.close();
					return false;
				}
				stmt = con.createStatement();
				String sql = "INSERT INTO " + table_name + " (username, password) VALUES(" + "'" + user_name.replace("'", "") +
								"'" + ", '" + password.replace("'", "") + "')";
				stmt.executeUpdate(sql);
				System.out.println("Inserted user login in given database...");
			}catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			}finally{
		      //finally block used to close resources
		      try{
		         if(stmt!=null)
		            con.close();
		      }catch(SQLException se){
		      }// do nothing
		      try{
		         if(con!=null)
		            con.close();
		      }catch(SQLException se){
		         se.printStackTrace();
		      }
			}
			return true;
		}
		return false;
	}
	public boolean confirmUser(String table_name, String username, String password){
		Connection con = null;
		Statement stmt = null;
		if(tableExists(table_name)){
			boolean check = false;
			try{
				con = DriverManager.getConnection(connectionUrl);
				String exists = "SELECT password FROM " + table_name + " WHERE username='" + username.replace("'", "") + "'";
				stmt = con.createStatement();
				ResultSet rs = stmt.executeQuery(exists);
				while(rs.next()){
					check = true;
					String confirm = rs.getString("password");
					if(!confirm.equals(password.replace("'", ""))){
						con.close();
						return false;
					}
				}
			}catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			}finally{
		      //finally block used to close resources
		      try{
		         if(stmt!=null)
		            con.close();
		      }catch(SQLException se){
		      }// do nothing
		      try{
		         if(con!=null)
		            con.close();
		      }catch(SQLException se){
		         se.printStackTrace();
		      }
			}
			return check;
		}
		return false;
	}
	public void createRestaurantTable(String table_name){
		Connection con = null;
		Statement stmt = null;
		try{
			con = DriverManager.getConnection(connectionUrl);
			stmt = con.createStatement();
			String sql = "CREATE TABLE " + table_name + " (id INTEGER not NULL, name VARCHAR(255), " +
							"address VARCHAR(255), number VARCHAR(255), categories VARCHAR(255), price INTEGER, " + 
							"times VARCHAR(255), url VARCHAR(255), icon VARCHAR(255), mapaddress VARCHAR(255))";
			stmt.executeUpdate(sql);
			System.out.println("Created table in given database...");
		}catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}finally{
		      //finally block used to close resources
		      try{
		         if(stmt!=null)
		            con.close();
		      }catch(SQLException se){
		      }// do nothing
		      try{
		         if(con!=null)
		            con.close();
		      }catch(SQLException se){
		         se.printStackTrace();
		      }
		}
	}
	public void insertRestaurants(String table_name,ArrayList<Restaurant> arr){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "INSERT INTO " + table_name + "(id, name, address, number, categories, price, times, url, icon, mapaddress) " +
							"VALUES(";
				for(int i = 0;i < arr.size();i++){
					Restaurant r = arr.get(i);
					String name = r.getName().replace("'", "");
					String address = r.getAddress().replace("'", "").replace("\"", "");
					String number = r.getNumber().replace("'", "");
					String categories = r.getCategories().replace("'", "").replace("\"", "");
					int price = r.getPrice();
					String times = r.getTimes().replace("'", "").replace("\"", "");
					String url = r.getUrl().replace("'", "").replace("\"", "");
					String icon = r.getIcon().replace("'", "");
					if(icon.equals("food")) icon = "foodpic";
					else if(icon.equals("dessert")) icon = "icecream";
					else if(icon.equals("drink")) icon = "hotdrink";
					String mapAddress = r.getMapAddress().replace("'", "");
					String values = Integer.toString(i) + ", " + "'" + name + "'" + 
									", " + "'" + address + "'" + 
									", " + "'" + number + "'" + 
									", " + "'" + categories + "'" +
									", " + Integer.toString(price) +
									", " + "'" + times + "'" +
									", " + "'" + url + "'" +
									", " + "'" + icon + "'" + 
									", " + "'" + mapAddress + "'" + 
									")";
					stmt.executeUpdate(sql + values);
				}
				System.out.println("Inserted restaurants in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public ArrayList<Restaurant> getRestaurantsFiltered(String table_name, ArrayList<String> categories, ArrayList<Integer> prices){
		Connection con = null;
		Statement stmt = null;
		ArrayList<Restaurant> ret = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name;
				if(categories.size() > 0){
					sql += " WHERE (categories LIKE '%";
					for(int i = 0;i < categories.size() - 1;i++){
						sql += categories.get(i) + "%" + "' OR categories LIKE '%";
					}
					sql += categories.get(categories.size() - 1) + "%')";
					if(prices.size() > 0){
						sql += " AND (price=";
						for(int i = 0;i < prices.size() - 1;i++){
							sql += Integer.toString(prices.get(i)) + " OR price=";
						}
						sql += Integer.toString(prices.get(prices.size() - 1)) + ")";
					}
				}
				else if(prices.size() > 0){
					sql += " WHERE (price=";
					for(int i = 0;i < prices.size() - 1;i++){
						sql += Integer.toString(prices.get(i)) + " OR price=";
					}
					sql += Integer.toString(prices.get(prices.size() - 1)) + ")";
				}
				ResultSet rs = stmt.executeQuery(sql);
				ret = new ArrayList<Restaurant>();
				while(rs.next()){
					Restaurant temp = new Restaurant();
					temp.setName(rs.getString(2));
					temp.setAddress(rs.getString(3));
					temp.setNumber(rs.getString(4));
					temp.setCategories(rs.getString(5));
					temp.setPrice(rs.getInt(6));
					temp.setTimes(rs.getString(7));
					temp.setUrl(rs.getString(8));
					temp.setIcon(rs.getString(9));
					temp.setMapAddress(rs.getString(10));
					ret.add(temp);
				}
				System.out.println("Selected restaurants in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return ret;
	}
	public ArrayList<Restaurant> getFourRestaurants(String table_name,ArrayList<String> categories,ArrayList<Integer> prices){
		ArrayList<Restaurant> ret = getRestaurantsFiltered(table_name, categories, prices);
		ArrayList<Double> means = getMeans("MEANS");
		if(ret.size() > 0){
			ArrayList<String> restaurant_names = new ArrayList<String>();
			for(Restaurant r : ret){
				restaurant_names.add(r.getName());
			}
			HashMap<Restaurant,Integer> hm = getRestaurantNumbers(table_name, restaurant_names);
			ArrayList<Integer> intarr = new ArrayList<Integer>(hm.values());
			ArrayList<Restaurant> restarr = new ArrayList<Restaurant>(hm.keySet());
			ArrayList<Double> dubarr = new ArrayList<Double>();
			for(int i = 0;i < intarr.size();i++){
				dubarr.add(means.get(intarr.get(i)));
			}
			quickSort(dubarr, restarr, 0, dubarr.size());
			ArrayList<Restaurant> topFour = new ArrayList<Restaurant>();
			if(restarr.size() > 3){
				for(int i = restarr.size() - 1;i >= restarr.size() - 4;i--){
					topFour.add(restarr.get(i));
				}
			}
			else{
				for(int i = restarr.size() - 1;i >= 0;i--){
					topFour.add(restarr.get(i));
				}
				ArrayList<Integer> indices = new ArrayList<Integer>();
				ArrayList<Double> doubles = new ArrayList<Double>();
				for(int i = 0;i < means.size();i++){
					if(indices.size() == 4 - restarr.size()){
						for(int j = 0;j < 4 - restarr.size();j++){
							if(Double.compare(means.get(i), doubles.get(j)) > 0){
								indices.add(j,i);
								doubles.add(j,means.get(i));
								indices.remove(4 - restarr.size());
								doubles.remove(4 - restarr.size());
								break;
							}
						}
					}
					else if(indices.size() > 0){
						boolean not_added = true;
						for(int j = 0;j < indices.size();j++){
							if(Double.compare(means.get(i), doubles.get(j)) > 0){
								indices.add(j,i);
								doubles.add(j,means.get(i));
								not_added = false;
							}
						}
						if(not_added) {
							indices.add(i);
							doubles.add(means.get(i));
						}
					}
					else{
						indices.add(i);
						doubles.add(means.get(i));
					}
				}
				ArrayList<Restaurant> remains = getRestaurantsByNumbers(table_name,indices);
				for(Restaurant r : remains){
					topFour.add(r);
				}
			}
			return topFour;
		}
		else{
			ArrayList<Integer> indices = new ArrayList<Integer>();
			ArrayList<Double> doubles = new ArrayList<Double>();
			for(int i = 0;i < means.size();i++){
				if(indices.size() == 4){
					for(int j = 0;j < 4;j++){
						if(Double.compare(means.get(i), doubles.get(j)) > 0){
							indices.add(j,i);
							doubles.add(j,means.get(i));
							indices.remove(4);
							doubles.remove(4);
							break;
						}
					}
				}
				else if(indices.size() > 0){
					boolean not_added = true;
					for(int j = 0;j < indices.size();j++){
						if(Double.compare(means.get(i), doubles.get(j)) > 0){
							indices.add(j,i);
							doubles.add(j,means.get(i));
							not_added = false;
						}
					}
					if(not_added) {
						indices.add(i);
						doubles.add(means.get(i));
					}
				}
				else{
					indices.add(i);
					doubles.add(means.get(i));
				}
			}
			ArrayList<Restaurant> topFour = getRestaurantsByNumbers(table_name,indices);
			return topFour;
		}
	}
	public ArrayList<Restaurant> getRestaurantsByNumbers(String table_name,ArrayList<Integer> indices){
		Connection con = null;
		Statement stmt = null;
		ArrayList<Restaurant> arr = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name + " WHERE (id='";
				for(int i = 0;i < indices.size() - 1;i++){
					sql += indices.get(i) + "' OR id='";
				}
				sql += indices.get(indices.size() - 1) + "')";
				ResultSet rs = stmt.executeQuery(sql);
				arr = new ArrayList<Restaurant>();
				while(rs.next()){
					Restaurant temp = new Restaurant();
					temp.setName(rs.getString(2));
					temp.setAddress(rs.getString(3));
					temp.setNumber(rs.getString(4));
					temp.setCategories(rs.getString(5));
					temp.setPrice(rs.getInt(6));
					temp.setTimes(rs.getString(7));
					temp.setUrl(rs.getString(8));
					temp.setIcon(rs.getString(9));
					temp.setMapAddress(rs.getString(10));
					arr.add(temp);
				}
				System.out.println("Selected user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return arr;
	}
	public HashMap<Restaurant,Integer> getRestaurantNumbers(String table_name,ArrayList<String> restaurants){
		Connection con = null;
		Statement stmt = null;
		HashMap<Restaurant,Integer> hm = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name + " WHERE (name='";
				for(int i = 0;i < restaurants.size() - 1;i++){
					sql += restaurants.get(i) + "' OR name='";
				}
				sql += restaurants.get(restaurants.size() - 1) + "')";
				ResultSet rs = stmt.executeQuery(sql);
				hm = new HashMap<Restaurant,Integer>();
				while(rs.next()){
					Restaurant temp = new Restaurant();
					temp.setName(rs.getString(2));
					temp.setAddress(rs.getString(3));
					temp.setNumber(rs.getString(4));
					temp.setCategories(rs.getString(5));
					temp.setPrice(rs.getInt(6));
					temp.setTimes(rs.getString(7));
					temp.setUrl(rs.getString(8));
					temp.setIcon(rs.getString(9));
					temp.setMapAddress(rs.getString(10));
					hm.put(temp, rs.getInt(1));
				}
				System.out.println("Selected user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return hm;
	}
	public void dropTable(String table_name){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "DROP TABLE " + table_name + " ";
				stmt.executeUpdate(sql);
				System.out.println("Deleted table in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public DoubleMatrix selectMatrix(String table_name, int num_Users) throws InterruptedException, ExecutionException{
		Connection con = null;
		Statement stmt = null;
		DoubleMatrix X = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name;
				ResultSet rs = stmt.executeQuery(sql);
				while(rs.next()){
					ResultSetMetaData rsmd = rs.getMetaData();
					int num_columns = rsmd.getColumnCount();
					DoubleMatrix temp = DoubleMatrix.zeros(1,num_columns - 2);
					if(X == null) X = DoubleMatrix.zeros(num_Users,num_columns - 2);
					int row = rs.getInt(1);
					for(int i = 3;i <= num_columns;i++){
						temp.put(0, i - 3, rs.getDouble(i));
					}
					X.putRow(row,temp);
				}
				int num_cols = X.getColumns();
				int[] cols = new int[num_cols];
				for(int i = 0;i < num_cols;i++){
					cols[i] = i;
				}
				List<HashMap<Double,DoubleMatrix>> lm = runAllMeans(X,cols);
				X = null;
				ArrayList<Double> arr = new ArrayList<Double>();
				for(int i = 0;i < lm.size();i++){
					HashMap<Double,DoubleMatrix> temp = lm.get(i);
					//assert one double key
					for(Double d : temp.keySet()){
						arr.add(d);
						if(X == null) X = temp.get(d);
						else X = DoubleMatrix.concatHorizontally(X, temp.get(d));
					}
				}
				String[] colstr = new String[arr.size()];
				for(int i = 0;i < arr.size();i++){
					colstr[i] = Integer.toString(i);
				}
				dropTable("MEANS");
				createTableMeans(colstr,"MEANS");
				insertTableMeans(arr,"MEANS");
				System.out.println("Selected in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			      
			}
		return X;
	}
	public int getUserNumber(String table_name, String user_name){
		Connection con = null;
		Statement stmt = null;
		int id = 0;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT id FROM " + table_name + " WHERE client=" + "'" + user_name.replace("'", "") + "'";
				ResultSet rs = stmt.executeQuery(sql);
				while(rs.next()){
					id = rs.getInt("id");
				}
				System.out.println("Selected user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return id;
	}
	public void updateUser(String table_name,String user_name,HashMap<Integer,Double> vals){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "UPDATE " + table_name + " SET ";
				String where = " WHERE client=" + "'" + user_name.replace("'","") +"'";
				for(Integer i : vals.keySet()){
					sql += "col" + Integer.toString(i) + "=" + Double.toString(vals.get(i)) + ", ";
				}
				sql = sql.substring(0,sql.length() - 2) + where;
				stmt.executeUpdate(sql);
				System.out.println("Updated user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	public ArrayList<Double> getMeans(String table_name){
		Connection con = null;
		Statement stmt = null;
		ArrayList<Double> ret = new ArrayList<Double>();
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name;
				ResultSet rs = stmt.executeQuery(sql);
				while(rs.next()){
					ResultSetMetaData rsmd = rs.getMetaData();
					int num_columns = rsmd.getColumnCount();
					for(int i = 1;i <= num_columns;i++){
						ret.add(rs.getDouble(i));
					}
				}
				System.out.println("Selected means in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return ret;
	}
	public ArrayList<Restaurant> getRestaurants(String table_name){
		Connection con = null;
		Statement stmt = null;
		ArrayList<Restaurant> arr = new ArrayList<Restaurant>();
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name;
				ResultSet rs = stmt.executeQuery(sql);
				while(rs.next()){
					Restaurant r = new Restaurant();
					r.setName(rs.getString(2));
					r.setAddress(rs.getString(3));
					r.setNumber(rs.getString(4));
					r.setCategories(rs.getString(5));
					r.setPrice(rs.getInt(6));
					r.setTimes(rs.getString(7));
					r.setUrl(rs.getString(8));
					r.setIcon(rs.getString(9));
					r.setMapAddress(rs.getString(10));
					arr.add(r);
				}
				System.out.println("Selected user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return arr;
	}
	public ArrayList<Double> getUser(String table_name, String user_name){
		Connection con = null;
		Statement stmt = null;
		ArrayList<Double> ret = new ArrayList<Double>();
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name + " WHERE client=" + "'" + user_name + "'";
				ResultSet rs = stmt.executeQuery(sql);
				while(rs.next()){
					ResultSetMetaData rsmd = rs.getMetaData();
					int num_columns = rsmd.getColumnCount();
					for(int i = 3;i <= num_columns;i++){
						ret.add(i - 3, rs.getDouble(i));
					}
				}
				ArrayList<Double> means = getMeans("MEANS");
				for(int i = 0;i < ret.size();i++){
					if(Double.compare(ret.get(i), 0.0) == 0){
						ret.set(i, means.get(i));
					}
				}
				System.out.println("Selected user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return ret;
	}
	public ArrayList<Double> getUserID(String table_name, int id){
		Connection con = null;
		Statement stmt = null;
		ArrayList<Double> ret = new ArrayList<Double>();
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "SELECT * FROM " + table_name + " WHERE id=" + "" + Integer.toString(id) + "";
				ResultSet rs = stmt.executeQuery(sql);
				while(rs.next()){
					ResultSetMetaData rsmd = rs.getMetaData();
					int num_columns = rsmd.getColumnCount();
					for(int i = 3;i <= num_columns;i++){
						ret.add(i - 3, rs.getDouble(i));
					}
				}
				System.out.println("Selected user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
		return ret;
	}
	public void deleteUser(String table_name, String user_name){
		Connection con = null;
		Statement stmt = null;
		try {
				con = DriverManager.getConnection(connectionUrl);
				stmt = con.createStatement();
				String sql = "DELETE FROM " + table_name + " WHERE client='" + user_name + "'";
				stmt.executeUpdate(sql);
				System.out.println("Deleted user in given database...");
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}finally{
			      //finally block used to close resources
			      try{
			         if(stmt!=null)
			            con.close();
			      }catch(SQLException se){
			      }// do nothing
			      try{
			         if(con!=null)
			            con.close();
			      }catch(SQLException se){
			         se.printStackTrace();
			      }
			}
	}
	private List<HashMap<Double,DoubleMatrix>> runAllMeans(DoubleMatrix X, int[] cols) 
			throws InterruptedException, ExecutionException {
		int threads = Runtime.getRuntime().availableProcessors();
		final int[] c = cols;
		final DoubleMatrix this_x = X;
		ExecutorService service = Executors.newFixedThreadPool(threads);
		List<Future<HashMap<Double,DoubleMatrix>>> futures = new ArrayList<Future<HashMap<Double,DoubleMatrix>>>();
		for (final int i : c) {
			Callable<HashMap<Double,DoubleMatrix>> callable = new Callable<HashMap<Double,DoubleMatrix>>() {
	            public HashMap<Double,DoubleMatrix> call() throws Exception {
	            	HashMap<Double,DoubleMatrix> output = runMean(this_x,i);
	                return output;
	            }
			};
			futures.add(service.submit(callable));
		}
		service.shutdown();
		List<HashMap<Double,DoubleMatrix>> outputs = new ArrayList<HashMap<Double,DoubleMatrix>>();
	    for (Future<HashMap<Double,DoubleMatrix>> future : futures) {
	        outputs.add(future.get());
	    }
	    return outputs;
	}
	private HashMap<Double,DoubleMatrix> runMean(DoubleMatrix X,int i){
		DoubleMatrix m = X.getColumn(i);
		int[] indices_nonzeros = m.findIndices();
		int num_rows = m.getRows();
		int num_nonzeros = indices_nonzeros.length;
		HashMap<Integer,Boolean> hm = new HashMap<Integer,Boolean>();
		for(int j = 0;j < num_nonzeros;j++){
			hm.put(indices_nonzeros[j], true);
		}
		int[] indices_zeros = new int[num_rows - num_nonzeros];
		int count = 0;
		for(int j = 0;j < num_rows;j++){
			if(!hm.containsKey(j)){
				indices_zeros[count] = j;
				count++;
			}
		}
		double mean = 0;
		if(indices_nonzeros.length > 0) mean = m.sum()/(indices_nonzeros.length);
		DoubleMatrix ret = m.put(indices_zeros, mean);
		HashMap<Double,DoubleMatrix> retmap = new HashMap<Double,DoubleMatrix>();
		retmap.put(mean, ret);
		return retmap;
	}/*
	private static HashMap<String,DoubleMatrix> meansMap(HashMap<String,ArrayList<Double>> hm) throws InterruptedException, ExecutionException{
		DoubleMatrix X = null;
		int count = 0;
		HashMap<Integer,String> stringmap = new HashMap<Integer,String>();
		for(String k : hm.keySet()){
			stringmap.put(count, k);
			count++;
			ArrayList<Double> temp = hm.get(k);
			int num_cols = temp.size();
			DoubleMatrix tm = DoubleMatrix.zeros(1,num_cols);
			for(int i = 0;i < num_cols;i++){
				tm.put(0, i, temp.get(i));
			}
			if(X == null) X = tm;
			else X = DoubleMatrix.concatVertically(X, tm);
		}
		int num_cols = X.getColumns();
		int[] cols = new int[num_cols];
		for(int i = 0;i < num_cols;i++){
			cols[i] = i;
		}
		List<DoubleMatrix> arr = runAllMeans(X,cols);
		X = arr.get(0);
		for(int i = 1;i < arr.size();i++){
			X = DoubleMatrix.concatHorizontally(X,arr.get(i));
		}
		HashMap<String,DoubleMatrix> ret = new HashMap<String,DoubleMatrix>();
		for(int i = 0;i < count;i++){
			String k = stringmap.get(i);
			ret.put(k, X.getRow(i));
		}
		return ret;
	}*/
}
