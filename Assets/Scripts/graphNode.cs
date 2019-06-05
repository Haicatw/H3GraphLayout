using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class graphNode : MonoBehaviour
{
    public int id;
    public Vector3 startMarker;
    public Vector3 endMarker;
    public float lerpTime;
    private Vector3 errorValue;
    private float lerpStartTime;

    // Start is called before the first frame update
    void Start()
    {
        startMarker = new Vector3(0, 0, 0);
        endMarker = new Vector3(-100, -100, -100);
        errorValue = new Vector3(-100, -100, -100);
    }

    // Update is called once per frame
    void Update()
    {
        if (transform.position != endMarker && endMarker != errorValue)
        {
            translateByanimation();
        }
        else
        {
            lerpStartTime = Time.time;
        }
    }

    public void setId(int id)
    {
        this.id = id;
    }

    public int getId()
    {
        return this.id;
    }

    public Vector3 lerp(float lerpStartTime)
    {
        float timeSinceStart = Time.time - lerpStartTime;
        float progressPercentage = timeSinceStart / lerpTime;
        return Vector3.Lerp(startMarker, endMarker, progressPercentage);
    }

    public void translateByanimation()
    {
        //float lerpStartTime = Time.time;
        transform.position = lerp(lerpStartTime);
    }

    public void setStartMarker(Vector3 startMarker)
    {
        this.startMarker = startMarker;
    }

    public void setEndMarker(Vector3 endMarker)
    {
        this.endMarker = endMarker;
    }

    public void setLerpTime(float lerpTime)
    {
        this.lerpTime = lerpTime;
    }
}
