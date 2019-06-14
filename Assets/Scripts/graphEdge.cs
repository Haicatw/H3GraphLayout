using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class graphEdge : MonoBehaviour
{
    public int parentNodeId;
    public int childNodeId;
    public GraphLayout graphLayoutScript;
    public bool initialized;
    private LineRenderer lr;

    // Start is called before the first frame update
    void Start()
    {
        GameObject graph = GameObject.FindGameObjectWithTag("GlobalManager");
        initialized = false;
        graphLayoutScript = graph.GetComponent<GraphLayout>();
    }

    // Update is called once per frame
    void Update()
    {
        if (initialized)
        {
            lr = this.GetComponent<LineRenderer>();
            lr.SetPosition(0, this.graphLayoutScript.nodesPrimitives[parentNodeId].transform.position);
            lr.SetPosition(1, this.graphLayoutScript.nodesPrimitives[childNodeId].transform.position);
        }
        
    }

    public void setParentId(int id)
    {
        this.parentNodeId = id;
    }

    public void setChildId(int id)
    {
        this.childNodeId = id;
    }

    public void setLineWidth(float lineWidth)
    {
        lr = this.GetComponent<LineRenderer>();
        lr.SetWidth(lineWidth, lineWidth);
    }

    public void setStartEndPosition(Vector3 start, Vector3 end)
    {
        lr.SetPosition(0, start);
        lr.SetPosition(1, end);
    }
}