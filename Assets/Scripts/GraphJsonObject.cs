using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class GraphJsonObject {
    public class JSONNode
    {
        public string id;
        public int level;
        public bool isRootNode;
    }

    public class JSONFullLink
    {
        public JSONNode source;
        public JSONNode target;
        public int value;
        public bool isTreeLink;
    }

    public class JSONTreeLink
    {
        public JSONNode source;
        public JSONNode target;
        public int value;
    }

    public class JSONMetadata
    {
        public int NumNodes;
        public int NumLinks;
        public int MaxLevel;
    }

    class JSONGraph
    {
        public List<JSONNode> jsonNodes;
        public List<JSONFullLink> jsonFullLinks;
        public List<JSONTreeLink> jsonTreeLinks;
        public JSONMetadata jsonMetadata;

        JSONGraph()
        {
            this.jsonNodes = new List<JSONNode>();
            this.jsonFullLinks = new List<JSONFullLink>();
            this.jsonTreeLinks = new List<JSONTreeLink>();
            this.jsonMetadata = new JSONMetadata();
        }
    }
}
