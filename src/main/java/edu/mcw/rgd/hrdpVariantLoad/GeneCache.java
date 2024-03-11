package edu.mcw.rgd.hrdpVariantLoad;

import javax.sql.DataSource;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 3/10/14
 * Time: 12:38 PM
 * <p>
 * gene rgd ids searchable by position and chromosome
 */
public class GeneCache {

    List<GeneCacheEntry> entries = new ArrayList<>();

    public int loadCache(int mapKey, String chromosome, DataSource ds) throws SQLException {

        entries.clear();

        String sql = "SELECT md.rgd_id,md.start_pos,md.stop_pos "+
	        "FROM maps_data md, rgd_ids r, genes g "+
			"WHERE md.chromosome=? AND md.map_key=? "+
			"AND md.rgd_id = r.RGD_ID and r.rgd_id = g.rgd_id "+
			"and r.OBJECT_STATUS = 'ACTIVE' and r.OBJECT_KEY = 1 "+
            "ORDER BY start_pos, stop_pos";
        Connection conn = ds.getConnection();
        PreparedStatement ps = conn.prepareStatement(sql);
        ps.setString(1, chromosome);
        ps.setInt(2, mapKey);
        ResultSet rs = ps.executeQuery();
        while( rs.next() ) {
            entries.add(new GeneCacheEntry(rs.getInt(1), rs.getInt(2), rs.getInt(3)));
        }
        conn.close();

        return entries.size();
    }

    /** return list of rgd ids for genes that match a snv variant given its position
     *
     * @param pos variant position
     * @return list of matching gene rgd ids, possibly empty
     */
    List<Integer> getGeneRgdIds(int start, int stop) {
        GeneCacheEntry key = new GeneCacheEntry(0, start, stop);
       /* int i = Collections.binarySearch(entries, key, new Comparator<GeneCacheEntry>() {
            public int compare(GeneCacheEntry o1, GeneCacheEntry o2) {
                System.out.println(o1.rgdId + "," + o2.rgdId);
                if( o2.stopPos < o1.startPos ) //5982114 < 5977943
                    return 1;
                if( o1.stopPos < o2.startPos ) //5979905 < 5982114
                    return -1;

                return 0;
            }
        });

        if( i<0 )
            return Collections.emptyList();
*/
        // there is a hit! add the hit to the results
        List<Integer> results = new ArrayList<Integer>();
        for(int i=0;i<entries.size();i++){
            GeneCacheEntry entry = entries.get(i);
            if(key.startPos >= entry.startPos && key.stopPos <= entry.stopPos)
                results.add(entry.rgdId);
            if (key.startPos <= entry.startPos && key.stopPos >= entry.startPos)
                results.add(entry.rgdId);
            if (key.startPos <= entry.stopPos && key.stopPos >= entry.stopPos)
                results.add(entry.rgdId);
        }
 /*       results.add(entries.get(i).rgdId);
        // look for possible other hits to the left of the hit index
        for( int j=i-1; j>=0; j-- ) {
            GeneCacheEntry e = entries.get(j);
            if( e.stopPos < pos )
                break;
            results.add(e.rgdId);
        }
        // look for possible other hits to the right of the hit index
        for( int j=i+1; j<entries.size(); j++ ) {
            GeneCacheEntry e = entries.get(j);
            if( e.startPos > pos )
                break;
            results.add(e.rgdId);
        }
        */
        return results;
    }

    class GeneCacheEntry {
        public int rgdId;
        public int startPos;
        public int stopPos;

        public GeneCacheEntry(int rgdId, int startPos, int stopPos) {
            this.rgdId = rgdId;
            this.startPos = startPos;
            this.stopPos = stopPos;
        }
    }
}
