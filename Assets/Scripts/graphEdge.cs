using System.Collections;
using System.Collections.Generic;
using UnityEngine;
/*
public class graphEdge : MonoBehaviour
{
    public Vector3 startMarkerStart;
    public Vector3 endMarkerStart;
    public Vector3 startMarkerEnd;
    public Vector3 endMarkerEnd;
    public int parentNodeId;
    public int childNodeId;
    public float lerpTime;
    private Vector3 errorValue;
    private float lerpStartTime;
    private LineRenderer lr;

    // Start is called before the first frame update
    void Start()
    {
        startMarkerStart = new Vector3(0, 0, 0);
        endMarkerStart = new Vector3(-100, -100, -100);
        startMarkerEnd = new Vector3(0, 0, 0);
        endMarkerEnd = new Vector3(-100, -100, -100);
        errorValue = new Vector3(-100, -100, -100);
        //lr = this.GetComponent<LineRenderer>();
    }

    // Update is called once per frame
    void Update()
    {
        lr = this.GetComponent<LineRenderer>();
        if (lr.GetPosition(0) != endMarkerStart 
            && endMarkerStart != errorValue 
            && lr.GetPosition(1) != endMarkerEnd
            && endMarkerEnd != errorValue)
        {
            translateByanimation();
        }
        else
        {
            lerpStartTime = Time.time;
        }
    }

    public Vector3 lerp(float lerpStartTime, Vector3 startMarker, Vector3 endMarker)
    {
        float timeSinceStart = Time.time - lerpStartTime;
        float progressPercentage = timeSinceStart / lerpTime;
        return Vector3.Lerp(startMarker, endMarker, progressPercentage);
    }

    public void translateByanimation()
    {
        //float lerpStartTime = Time.time;
        lr = this.GetComponent<LineRenderer>();
        lr.SetPosition(0, lerp(lerpStartTime, startMarkerStart, endMarkerStart));
        lr.SetPosition(1, lerp(lerpStartTime, startMarkerEnd, endMarkerEnd));
    }

    public void setStartMarkerStart(Vector3 startMarker)
    {
        this.startMarkerStart = startMarker;
    }

    public void setEndMarkerStart(Vector3 endMarker)
    {
        this.endMarkerStart = endMarker;
    }

    public void setStartMarkerEnd(Vector3 startMarker)
    {
        this.startMarkerEnd = startMarker;
    }

    public void setEndMarkerEnd(Vector3 endMarker)
    {
        this.endMarkerEnd = endMarker;
    }

    public void setLerpTime(float lerpTime)
    {
        this.lerpTime = lerpTime;
    }

    public void setMaterial(Material mat, Color color)
    {
        lr.material = mat;
        lr.SetColors(color, color);
    }

    public void setLineWidth(float lineWidth)
    {
        lr = this.GetComponent<LineRenderer>();
        lr.SetWidth(lineWidth, lineWidth);
    }

    public void setParentId(int id)
    {
        this.parentNodeId = id;
    }

    public void setChildId(int id)
    {
        this.childNodeId = id;
    }

    public void setStartEndPosition(Vector3 start, Vector3 end)
    {
        lr.SetPosition(0, start);
        lr.SetPosition(1, end);
    }
}
*/
public class graphEdge : MonoBehaviour
{
    public int parentNodeId;
    public int childNodeId;
    public GraphLayout graphLayoutScript;
    private LineRenderer lr;

    // Start is called before the first frame update
    void Start()
    {
        GameObject graph = GameObject.FindGameObjectWithTag("GlobalManager");
        graphLayoutScript = graph.GetComponent<GraphLayout>();
    }

    // Update is called once per frame
    void Update()
    {
        lr = this.GetComponent<LineRenderer>();
        lr.SetPosition(0, this.graphLayoutScript.nodesPrimitives[parentNodeId].transform.position);
        lr.SetPosition(1, this.graphLayoutScript.nodesPrimitives[childNodeId].transform.position);
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
}